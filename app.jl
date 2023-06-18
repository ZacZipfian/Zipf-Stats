module App
# set up Genie development environmet
using GenieFramework
using StatsBase, DataFrames, Distributions
@genietools

Genie.config.cors_headers["Access-Control-Allow-Origin"]  =  "*"
Genie.config.cors_headers["Access-Control-Allow-Headers"] = "Content-Type"
Genie.config.cors_headers["Access-Control-Allow-Methods"] = "GET,POST,PUT,DELETE,OPTIONS"
Genie.config.cors_allowed_origins = ["*"]

const FILE_PATH = "upload"
mkpath(FILE_PATH)

function TransMatxx(Markov)
    n = size(Markov,1);
    trans = zeros(n,n);
    for i = 1:n
       trans[i,:] = (Markov[i,:]*(1/sum(Markov[i,:])))    
    end
    return trans
end

function mc_path(P, char_set; init , sample_size)
    @assert size(P)[1] == size(P)[2] # square required
    N = size(P)[1] # should be square

    # create vector of discrete RVs for each row
    dists = [Categorical(P[i, :]) for i in 1:N]

    # setup the simulation
    X = fill(0, sample_size) # allocate memory, or zeros(Int64, sample_size)
    X[1] = init # set the initial state
    Random_text = string(char_set[init])

    for t in 2:sample_size
        dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
        X[t] = rand(dist) # draw new value
        Random_text = string(Random_text,"",char_set[X[t]])
    end
    return Random_text
end

function doc_processing(FILE_PATH, Document_name)
    #words
    DocuString = read(joinpath(FILE_PATH, Document_name), String)
    WordsString = replace(DocuString, r"[^A-Za-z]" => " ")
    WordsString = lowercase(WordsString)
    WordVector = split(WordsString, " ")  
    WordVector = filter!(e->e!="",WordVector)
    Word_count = length(WordVector)  
    Col_words = countmap(WordVector)
    unq_wrds = length(Col_words)     
    Col_words = sort(collect(Col_words), by=x->x[2], rev=true)
    words = [first(p) for p in Col_words] 
    countW = [last(p) for p in Col_words]
    distrW = [(last(p)/Word_count) for p in Col_words]
    df1 = DataFrame(row = log.(1:unq_wrds), word = words, count = log.(countW), distribution = distrW)

    #symbols
    DocuString = replace(DocuString, "\ufeff" => "")
    DocuString = replace(DocuString, "\r" => "")
    DocuString = replace(DocuString, "\n" => "") 
    DocuString = replace(DocuString, "\$" => "")
    symbols = split(DocuString, "")
    CSym = countmap(symbols)
    CSym = sort(collect(CSym), by=x->x[2], rev=true)
    symsl = [first(p) for p in CSym]
    countS = [last(p) for p in CSym] 
    S_dist = [(last(p)/(sum(countS))*100) for p in CSym]  
    df2 = DataFrame(row = 1:length(countS), symbol = symsl, count = countS, distribution = S_dist)
 
    #Markov chain
    plcS = [first(p) for p in CSym]
    n = length(countS)
    ValidIndx = collect(eachindex(DocuString))
    CPrev = DocuString[1]
    markov_chain = zeros(n,n)

    for i in 2:length(ValidIndx)
        CCur = DocuString[ValidIndx[i]]
        Index1 = findall( x -> x == string(CPrev), plcS )
        Index2 = findall( x -> x == string(CCur), plcS )
        markov_chain[Index1[1], Index2[1]] += 1
        CPrev = CCur
    end
    transition_matrix = TransMatxx(markov_chain)
    df3 = DataFrame(rotr90(reverse(transition_matrix, dims = 1)), :auto)
    df3 = rename!(df3, ["x$i" => namess for (i, namess) in enumerate(plcS)])
    df3[!, :symb] = plcS
    df3[!, :rank] = 1:(size(df2,1))

    #miscellaneous
    trw = 0
    for j in 1:unq_wrds
        trw = trw + (length(words[j]) * countW[j])
    end
    av_wrd_count = (trw / Word_count) 
    W_dist = [(last(p)/(Word_count) * log2(last(p)/(Word_count))) for p in Col_words]  
    shan_entropy = sum(W_dist)*-1 

    return Word_count, unq_wrds, av_wrd_count, shan_entropy, df1, df2, df3, transition_matrix
end

const plot_colors = ["#0779e4", "#525254", "#da7c2e"]

# add reactive code to make the UI interactive
@app begin
    # reactive variables are tagged with @in and @out
    @out upfiles = readdir(FILE_PATH)
    @in Document_sel = "hamlet_full.txt"
    @out Markov_opts = ["a"]
    @in markov_sel = "a"
    @out genText = " "
    @in markov_length = 20
    @out zipfplot = PlotData[]
    @out layzipf = PlotLayout( xaxis=[PlotLayoutAxis(xy="x", title = "Log(Rank)")], 
                               yaxis=[PlotLayoutAxis(xy="y", title = "Log(Frequency)")])
    @out symb_distr = PlotData[]
    @out laysymb = PlotLayout( xaxis=[PlotLayoutAxis(xy="x", title ="Rank")], 
                               yaxis=[PlotLayoutAxis(xy="y", title ="Percentage of Text")])
    @out markovprobdist = PlotData[]
    @out laymarkov = PlotLayout( xaxis=[PlotLayoutAxis(xy="x", title ="Rank")], 
                                 yaxis=[PlotLayoutAxis(xy="y", title ="Probability")])
    @out graphname = "The probability that a is followed by:"
    @in generateM = false
    @in clear = false
    @out msg = "Word count: 0, Unique words: 0, Average word length: 0, Text entropy: 0"
    
    # @private defines a non-reactive variable
    @private Word_count = 0.0
    @private unq_wrds = 0.0
    @private av_wrd_count = 0.0
    @private shan_entropy = 0.0
    @private allowgen = false

    # watch a variable and execute a block of code when its value changes
    @event :uploaded begin
        upfiles = readdir(FILE_PATH)
    end
    @onchange Document_sel begin
        upfiles = readdir(FILE_PATH)
        global Word_count, unq_wrds, av_wrd_count, shan_entropy, df1, df2, df3, transition_matrix = doc_processing(FILE_PATH, Document_sel)
        msg = "Word count: $Word_count, Unique words: $unq_wrds, Average word length: $av_wrd_count, Text entropy: $shan_entropy"
        zipfplot = plotdata(df1, :row, :count, marker = PlotDataMarker(color = plot_colors[1]), xaxis="x", yaxis="y"; groupfeature = :word)
        symb_distr = plotdata(df2, :row, :distribution, plot=StipplePlotly.Charts.PLOT_TYPE_BAR, marker = PlotDataMarker(color = plot_colors[2]); groupfeature = :symbol)
        Markov_opts = df2[!, :symbol]
        markovprobdist = plotdata(df3, :rank, Symbol(markov_sel), plot=StipplePlotly.Charts.PLOT_TYPE_BAR, marker = PlotDataMarker(color = plot_colors[3]); groupfeature = :symb)
        graphname = "The probability that $markov_sel is followed by:"
        allowgen = true
    end
    @onchange markov_sel begin
        markovprobdist = plotdata(df3, :rank, Symbol(markov_sel), plot=StipplePlotly.Charts.PLOT_TYPE_BAR, marker = PlotDataMarker(color = plot_colors[3]); groupfeature = :symb)
        graphname = "The probability that $markov_sel is followed by:"
    end
    @onchange generateM begin
        if allowgen == true
            genText = mc_path(transition_matrix, Markov_opts, init=rand(1:(size(df2,1))), sample_size=markov_length)
        end
    end
    @onchange clear begin
        n = (length(upfiles))
        for i in 1:n 
            del_file = upfiles[i]
            rm(joinpath(FILE_PATH, del_file))
        end
        upfiles = readdir(FILE_PATH)
    end

    route("/upload", method = POST) do
        files = Genie.Requests.filespayload()
        for f in files
            write(joinpath(FILE_PATH, f[2].name), f[2].data)
        end
        if length(files) == 0
            @info "No file uploaded"
        end
         upfiles = readdir(FILE_PATH)
        return "Upload finished"
    end
end

# loaded on access
@page("/", "app.jl.html")
Server.isrunning() || Server.up()
end
