#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 9/13/2015
   e-mail: tim.sterne.weiler@utoronto.ca
=#

include("dictext.jl")

using ArgParse
using DataFrames
using Distributions

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--bed"
      help = "toggle output as bed format (default on)"
      action = :store_true
    "--biot"
      help = "only output this specific biotype (default off)"
      arg_type = ASCIIString
    "--background"
      help = "gene expression levels background file {geneid first column, then data replicates cols}"
      arg_type = ASCIIString
    "--backheader"
      help = "is there a header in the background file? (default true)"
      action = :store_false
  end
  return parse_args(s)
end


###################################################################
function splitgb( coord )
  spl = split( coord, ':' )
  @assert( length(spl) > 1, "$coord does not split : into two!" )
  nums = map( x->parse(Int, x), split( spl[2], '-' ) )
  length(spl) >= 3 ? (spl[1], nums..., spl[3]) : (spl[1], nums...)
end


function read_lig_parse_reg( io; biot=nothing )
   const geneind = 25
   const biotind = 24
   const seqind = 11
   const posind = 17
   const refind = 14

   expfile = []
   bygenesym = Dict{ASCIIString,Int}()
   genesum = 0

   for l::ASCIIString in eachline( io )
      l[1] == 'S' || continue
      s = split(strip(l), '\t')
      g1,g2 = split( s[geneind], ':' )
      b1,b2 = split( s[biotind], ':' )
      #println("$g1 -- $g2, $b1 -- $b2")
      g1 == g2 || continue
      (biot != nothing && (b1 != biot || b2 != biot)) && continue
      s[posind][1:2] == "NA" && continue
      c1,c2 = map( x->splitgb(x), split( s[posind], ',' ) )
      ci1,ci2 = c1[2],c2[2]
      c1[1] == c2[1] || continue
      @assert( c1[end] == c2[end], "$(c1[end]) does not equal $(c2[end]) !!" )
      s1,s2 = replace(s[seqind], r"[a-z]", "") |> x->split( x, '_' )
      #println("$(s[refind]) is the ref ind")
      r1,r2 = split( s[refind], ',' ) |> x->tuple(map(y->parse(Int,y), x)...)
#=      if c1[end] == "+" # if positive
         ci1 += length(s1) - 1
      else # negative
         ci1 -= length(s1) + 1
      end=#
      first,second = tuple(sort([ci1,ci2])...)
      dist = abs( r1 - r2 )
      push!(expfile, (c1[1], first, second, g1, dist, c1[end]) )
      dinc!( bygenesym, ASCIIString(g1) )
      genesum += 1
   end
   expfile, bygenesym, genesum
end

function print_bed( io, reg )
   for bedtup in reg
      bedstr = join(bedtup, '\t')
      println( io , bedstr )
   end
end

function print_pvals( io, pvals::Dict{ASCIIString,Float64} )
   for k in keys(pvals)
      println( io, "$k\t$(pvals[k])" )
   end
end

function load_gene_exp( filename; header=true, ind=2 )
   df = readtable( filename, separator='\t', header=header)
   num_inds = get_numeric_cols( df )
   divsum( arr ) = arr / sum(arr)
   for i in num_inds
      df[i] = divsum( df[i] )
   end
   da = Vector{Float64}()
   for r in eachrow( df[n] ) 
      push!(da, mean(convert(Array, r)))
   end
   df[:mean] = divsum( da )
   geneDic = Dict{ASCIIString,Float64}()
   for i in 1:nrow(df)
      geneDic[ df[i,ind] ] = df[i,:mean]
   end
   geneDic
end

function get_numeric_cols( df::DataFrame )
   num_inds = Int[]
   for i in 1:length(df)
      if typeof(df[i][1]) <: Real
         push!(num_inds, i)
      end
   end
   num_inds
end

function main()
   pargs = parse_cmd()

   reg,genecnt,genesum = read_lig_parse_reg( STDIN, biot=pargs["biot"] )
   if pargs["bed"]
      #print_bed( STDOUT, reg )
   end

   if pargs["background"] != nothing
      back = load_gene_exp( pargs["background"], header=pargs["backheader"] )
      pvals = Dict{ASCIIString,Float64}()
      for k in keys(genecnt)
         haskey(back, k) || continue
         bin = Binomial(Int(genesum), back[k])
         pvals[ k ] = ccdf( bin, genecnt[k]-1 )
      end
      n = length(keys(pvals))
      for k in keys(pvals)
         pvals[k] = min( pvals[k] * n, 1 )
      end
      open( "n_vs_ge-back.pval", "r" ) do fh
         print_pvals( fh, pvals )
      end
   end   

end

main()
