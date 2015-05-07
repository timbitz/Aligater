#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 3/16/2015
   e-mail: tim.sterne.weiler@utoronto.ca
 
   dev on stable julia 0.3.8
=# 

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ## INITIALIZATION  ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

using ArgParse
using StatsBase
using Distributions

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--geneid", "-g"
      help = "flag to collapse geneids"
      arg_type = Bool
      default = false
    "--biotype", "-b"
      help = "flag to collapse biotypes"
      arg_type = Bool
      default = false
  end
  return parse_args(s)
end

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ###  FUNCTIONS ## ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #

include("dictext.jl") #--> dush!, dinc!, dnorm!

function isUniqueJunc!( used::Dict{ASCIIString,Bool}, seq::ASCIIString, genes )
  m = match(r"([AGCTUN]{5}_[AGCTUN]{5})", seq)
  cap = m.captures[1]
  m = match(r"([AGCTUN]+)_([AGCTUN]+)", seq)
  lLen,rLen = length(m.captures[1]), length(m.captures[2])
  geneid = sort( genes )
  key = join(geneid, ":") * "_$cap:$lLen:$rLen"
  retbool = false
  if !( haskey( used, key ) )
      used[key] = true
      retbool = true
  end
  retbool
end #--> Bool

function setClass!( doubleDict::Dict{ASCIIString,Dict{ASCIIString,Char}}, class::Char, seqs, genes; size=32 )
  @assert( length(seqs) == 2 && length(genes) == 2 )
  classHeirPair( a::Char, b::Char ) = a < b ? b : a  #--> Char
  classHeirArray( arr::Array{Char,1} ) = shift!( reverse( sort( arr ) ) ) #--> Char
  retarray = Char[]
  lenA,lenB = map(length, seqs)
  # pad a sequence if it is shorter than our target window size
  seqA = seqs[1] * repeat(".", max( size - lenA, 0 ))
  seqB = seqs[2] * repeat(".", max( size - lenB, 0 ))
  # reset the lengths after padding
  lenA,lenB = length(seqA),length(seqB)
  it=max(1, floor(size/2)) # calculate the step size, as half the window size
  # iterate through the cartesian product of windows in seqA and seqB by it step size
  for i = 1:it:(lenA-size)+1, j = 1:it:(lenB-size)+1
    winA = seqA[i:i+size-1] # access substrings 
    winB = seqB[j:j+size-1]
    keyA,keyB = sort( [winA,winB] )
    # use non-exported ht_keyindex function to avoid over indexing hash
    indxA = Base.ht_keyindex( doubleDict, keyA )
    if indxA <= 0 #initialize hash if necessary
      doubleDict[keyA] = Dict{ASCIIString,Char}()
      indxA = Base.ht_keyindex( doubleDict, keyA )
    end
    dictA = doubleDict.vals[indxA]
    indxB = Base.ht_keyindex( dictA, keyB )
    if indxB <= 0
      dictA[keyB] = class #set val if none exists
    else
      curVal = classHeirarchy( dictA.vals[indxB], class )
      dictA.vals[indxB] = curVal
      push!(retarray, curVal)
    end
  end
  classHeirArray( retarray )
end #--> Char

## ### ### ### ### ### ### ### ### ### ### ### ### ## #
# ### ### ### ### ### MAIN  ### ### ### ### ### ### ##
## ### ### ### ### ### ### ### ### ### ### ### ### ## #
function main()
  pargs = parse_cmd()

  # simplify ridiculous split string type
  # in case julia changes by version
  const stype = typeof( split("a:b", ':') )

  djunc = Dict{ASCIIString,Bool}()
  dstore = Dict{Char,Array{stype,1}}()

  seqInd = 11
  geneInd = 25

  for i in eachline( STDIN )

    s = split(chomp(i), '\t')
    genes = split(s[geneInd], ':')

    # test if this is a unique junction/readset
    if isUniqueJunc!( djunc, s[seqInd], genes )

      #set data
      @assert( length(s[1]) == 1 )
      curclass = Char( s[1] )
      dush!( dstore, curclass, s, arrType=typeof(s) )
      
    end
  end
end

main()
