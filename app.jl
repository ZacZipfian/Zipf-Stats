module App
# set up Genie development environmet
using GenieFramework
using StatsBase, DataFrames
@genietools

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

function doc_processing(FILE_PATH, Document_name)
    #words
    DocuString = read(joinpath(FILE_PATH, Document_name), String)
    WordsString = replace(DocuString, r"[^A-Za-z]" => " ")
    WordsString = lowercase(WordsString)
    WordVector = split(WordsString, " ")  
    WordVector = filter!(e->e!="",WordVector)
    Word_count = length(WordVector)  #return me 
    Col_words = countmap(WordVector)
    unq_wrds = length(Col_words)     #return me 
    Col_words = sort(collect(Col_words), by=x->x[2], rev=true)
    words = [first(p) for p in Col_words] #return me 
    countW = [last(p) for p in Col_words] #return me 
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
    symsl = [first(p) for p in CSym] #return me 
    countS = [last(p) for p in CSym] 
    S_dist = [last(p)/(sum(countS)) for p in CSym] #return me  
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
    transition_matrix = TransMatxx(markov_chain) #return me 

    #miscellaneous
    trw = 0
    for j in 1:unq_wrds
        trw = trw + (length(words[j]) * countW[j])
    end
    av_wrd_count = (trw / Word_count) #return me 

    W_dist = [(last(p)/(Word_count) * log2(last(p)/(Word_count))) for p in Col_words]  
    shan_entropy = sum(W_dist)*-1 #return me 
 
    return Word_count, unq_wrds, av_wrd_count, shan_entropy, df1, df2, transition_matrix

end

# add reactive code to make the UI interactive
@app begin
    Document_name = "hamlet_full.txt"
    @out Markov_opts = ["a", "b", "c"]
    @in markov_sel = ["a"]
    @out genText = " "
    @in markov_length = 0
    @out zipfplot = PlotData[]
    @out symb_distr = PlotData[]
    @out markovprobdist = PlotData[]


    # reactive variables are tagged with @in and @out
    @in generateM = false
    @in start = false
    @out msg = "Word count: 0, Unique words: 0, Average word length: 0, Text entropy: 0"
    # @private defines a non-reactive variable
    @private Word_count = 0.0
    @private unq_wrds = 0.0
    @private av_wrd_count = 0.0
    @private shan_entropy = 0.0

    # watch a variable and execute a block of code when
    # its value changes
    @onchange start begin
        Word_count, unq_wrds, av_wrd_count, shan_entropy, df1, df2, transition_matrix = doc_processing(FILE_PATH, Document_name)
        msg = "Word count: $Word_count, Unique words: $unq_wrds, Average word length: $av_wrd_count, Text entropy: $shan_entropy"
        zipfplot = plotdata(df1, :row, :count; groupfeature = :word)
        symb_distr = plotdata(df2, :row, :distribution, plot=StipplePlotly.Charts.PLOT_TYPE_BAR; groupfeature = :symbol)
        Markov_opts = df2[!, :symbol]
        indexii = findall( x -> x == string(markov_sel), Markov_opts )
        markovprobdist = plotdata(transition_matrix[indexii, :])
    end
end

# register a new route and the page that will begin
# loaded on access
@page("/", "app.jl.html")
Server.isrunning() || Server.up()
end
