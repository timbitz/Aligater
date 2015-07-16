#!/usr/bin/env julia
#=
   Author: Tim Sterne-Weiler, 6/23/2015
   e-mail: tim.sterne.weiler@utoronto.ca
=#

using ArgParse

function parse_cmd()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--snofile"
      help = "snoRNA fasta file"
      arg_type = ASCIIString
  end
  return parse_args(s)
end

##################################################################
function lenientregex( s::ASCIIString )
   levels = ASCIIString[]
   nlevels = 6
   for i in 1:nlevels
      push!(levels, s)  # initialize levels
   end 
   # this loop requires the first and last bases to match exactly
   for c in 2:(length(s)-1)
      pre = s[1:(c-1)]
      post = s[(c+1):end]
      levels[2] *= "|" * pre * "[atgcuATGCU]" * post
      levels[4] *= "|" * pre * "." * string(s[c]) * post
      levels[6] *= "|" * pre * "." * "[atgcuATGCU]" * post
   end
   # this one allows first and last bases to be lenient
   for c in 1:length(s)
      pre = s[1:(c-1)]
      post = s[(c+1):end]
      levels[3] *= "|" * pre * "[atgcuATGCU]" * post
      levels[5] *= "|" * pre * "." * string(s[c]) * post
   end   
   ret = Regex[]
   for i in 1:nlevels
     push!(ret, Regex(levels[i]))
   end
   ret
end #-->Array{Regex,1}

function searchregex( regexarray::Array{Regex,1}, s::ASCIIString )
   for r in regexarray
      m = match(r, s)
      if m != nothing
        return (m.offset,m.match)
      end
   end
   (0,"")
end #-->Tuple{Int64,ASCIIString}

function annotate_cdbox( sno::ASCIIString )
   cboxreg = lenientregex("TGATGA")
   dboxreg = lenientregex("CTGA")

   cboxseg = sno[1:20]
   cbox = searchregex(cboxreg, cboxseg)
   @assert(cbox[1] > 0, "Cannot find cbox in $cboxseg!")

   dboxseg = sno[(end-20):end]
   (dboxpos,dboxmatch) = searchregex(dboxreg, dboxseg)
   @assert(dboxpos > 0, "Cannot find dbox in $dboxseg!")
   dboxpos += length(sno) - 20
   dbox = (dboxpos, dboxmatch)

   dprimeseg = sno[(cbox[1]+16):(dboxpos-14)]
   (dprimepos, dprimematch) = searchregex(dboxreg, dprimeseg)
   dprimepos += cbox[1] + 16   
   dprime = (dprimepos, dprimematch)

   println(STDERR, "$(cbox[2]) ... $(dprime[2]) ... $(dbox[2])")

   (cbox, dprime, dbox)
end

function main()   
   for l::ASCIIString in eachline(STDIN)
      sno = uppercase(chomp(l))
            
   end
end
###################################################################
main()
