#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 9/13/2015
   e-mail: tim.sterne.weiler@utoronto.ca
=#

using ArgParse

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--bed"
      help = "toggle output as bed format (default on)"
      action = :store_true
    "--biot"
      help = "only output this specific biotype (default off)"
      arg_type = ASCIIString
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
      println("$(c1[1])\t$first\t$second\t$g1\t$dist\t$(c1[end])\t$s1\t$s2")
   end
end

function print_bed( io, reg )
   for bedtup in reg
      bedstr = join(bedtup, '\t')
      println( io , bedstr )
   end
end

function main()
   pargs = parse_cmd()

   reg = read_lig_parse_reg( STDIN, biot=pargs["biot"] )
   if pargs["bed"]
      #print_bed( STDOUT, reg )
   end
end

main()
