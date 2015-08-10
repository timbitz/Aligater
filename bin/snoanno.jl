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
    "--regex"
      help = "snoRNA fasta header regex with capture [optional]"
      arg_type = ASCIIString
      default = "^>\\s*.*_.*_(\\S+)_.*"
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
   #@assert(cbox[1] > 0, "Cannot find cbox in $cboxseg!")

   dboxseg = sno[(end-20):end]
   (dboxpos,dboxmatch) = searchregex(dboxreg, dboxseg)
   #@assert(dboxpos > 0, "Cannot find dbox in $dboxseg!")
   dboxpos += length(sno) - 20
   dbox = (dboxpos, dboxmatch)

   dprimeseg = sno[(cbox[1]+16):(dboxpos-14)]
   (dprimepos, dprimematch) = searchregex(dboxreg, dprimeseg)
   dprimepos += cbox[1] + 16   
   dprime = (dprimepos, dprimematch)

   println(STDERR, "$(cbox[2]) ... $(dprime[2]) ... $(dbox[2])")

   (cbox, dprime, dbox)
end #-->Tuple{Tuple}

# this function reads a fasta file/io and returns a dict of {name,seq}
function readfasta( io; regex = r">\s*(\S+)" )
   rethash = Dict{ASCIIString,ASCIIString}()  

   #internal funct for fasta header
   function checkname( head::ASCIIString, namereg::Regex )
      res = match(namereg, head)
      @assert(res != nothing, "$head looks to be an incorrectly formated fasta header!")
      chomp(res.captures[1])
   end
  
   curseq = ""
   head = readline(io)
   curname = checkname( head, regex ) 
   for line::ASCIIString in eachline( io )
      #finish up
      if line[1] == '>' #header line
         # push to dict
         rethash[curname] = curseq
         curname = checkname( line, regex )
         curseq = ""
      else #sequence line
        curseq *= chomp( uppercase(line) )
      end
   end
   rethash
end #--> Dict{ASCIIString,ASCIIString}

# find the distance from the 5' end of the transcript of the antisense region
function bind_distance( ind::Int64, ligstruct::String, startpos::Int64 )
   antiregex = r"([ATGCU]+(?:[\.\(\)]{1,5}[ATGCU]+)?)"
   m = match( antiregex, ligstruct ) # match antisense site
   cut = match( r"(\|)", ligstruct )
   m == nothing && return (1,"",1)
   cap = m.captures[1]
   len = length(cap)
   offset,lig = 0,0
   if ind <= 1
      offset = startpos + m.offsets[1]
      lig = startpos + (cut == nothing ? length(ligstruct) : cut.offsets[1])
   else
      barpos = cut == nothing ? 0 : cut.offsets[1]
      offset = startpos - barpos + m.offsets[1]
      lig = startpos
   end
   (offset, cap, lig)
end #--> Tuple{Int64,Int64,Int64}

function print_heatrow{I <: Integer}( io, rowname::ASCIIString, offset::I, antiseq::String, ligpos::I, anchorpos::I, totallen::I, boxes::Tuple; nais="NA", order="" )
   println(STDERR, "$io\t$rowname\t$offset\t$antiseq\t$ligpos\t$anchorpos")
   prestr = rowname * "\t"
   str = repeat( "0", totallen )
   for i = 1:3
     (off, seq) = boxes[i]
     off > 0 || continue
     asnum = replace( seq, r"[AUCGT]", string(i+2) )
     str = str[1:off-1] * asnum * str[off+length(asnum):end]
   end

   antinum = replace( antiseq, r"[AUCGT]", "1" ) |> x->replace( x, r"\.", "0" )
   str = str[1:offset-1] * antinum * str[offset+length(antinum):end]

   str = ligpos > totallen ? str : str[1:ligpos-1] * "2" * str[ligpos+1:end]
   rev = reverse(str)
   pad = totallen - anchorpos
   outiter = [ pad:-1:0; -1:-1:(-1*anchorpos)] 
   for i = 1:length(rev)
      cur = rev[i]
      println( io, prestr * string(outiter[i]) * order * "\t$cur" )
   end
end

###################################################################
function main()
   pargs = parse_cmd()
   
   fastareg = pargs["regex"] == nothing ? r"^>\s*(\S+)" : Regex(pargs["regex"])
   
   @assert(pargs["snofile"] != nothing, "you have to provide --snofile")
   snofile = pargs["snofile"]
   fh = open( snofile, "r" )
   snodict = readfasta( fh, regex=fastareg )
   close(fh)
   
   # try to annotate the C/D boxes
   cdboxhash = Dict{ASCIIString, Tuple}()
   for k in keys(snodict)
      cdboxhash[k] = annotate_cdbox( snodict[k] )
   end
  
   const genInd = 25

   output = Dict{ASCIIString, Tuple}()
   ligdist = Dict{ASCIIString, Int64}()
   capsize = Dict{ASCIIString, Int64}()
   annobox = Dict{ASCIIString, Tuple}()

   # now lets go through the pvl/lig data
   for l in eachline(STDIN)
      s = split( l, '\t' )
      @assert( ismatch(r"SNORD", s[genInd]), "$(s[genInd]) doesn't appear to have a snoRNA in it? (SNORD)\n" )
      ids = map(x->chomp(x), split( s[genInd], ':' ))
      ind = ismatch(r"SNORD", ids[1]) ? 1 : 2
      @assert( ismatch(r"SNORD", ids[ind]) )
      struct = s[ind+18]
      start = parse(Int, split( s[14], ',' )[ind])
      (antioffset, cap, lig) = bind_distance( ind, struct, start )
      @assert( haskey(cdboxhash, "SNORD9" ), "$(ids[ind]) doesn't have a key!" )
      cdboxes = cdboxhash[ ids[ind] ]
      snolen  = length( snodict[ ids[ind] ] )
      curname_key = ids[ind] * "_" * ids[ ind == 1 ? 2 : 1 ]
      if antioffset <= cdboxes[2][1] # D' box offset
         #curname_key = ids[ind] * "_DpBOX"
         outtup = (curname_key, antioffset, cap, lig, cdboxes[2][1], snolen ) 
         outdist = ( cdboxes[2][1] - (round(Int,length(cap) / 2) + antioffset) )
      elseif antioffset <= cdboxes[3][1] # D box offset
         #curname_key = ids[ind] * "_DBOX"
         outtup = (curname_key, antioffset, cap, lig, cdboxes[3][1], snolen ) 
         outdist = ( cdboxes[3][1] - (round(Int,length(cap) / 2) + antioffset) )
      else
         continue
      end

      # add to output dicts
      storedlen = get( capsize, curname_key, 0 )
      if length(cap) > storedlen
         output[curname_key]  = outtup # store anno tuple
         ligdist[curname_key] = outdist # store distance to nearest D box
         capsize[curname_key] = length(cap) # store length of antisense
         annobox[curname_key] = cdboxes #store tuple of c/d box annotations
      end
   end

   # print in order
   keysarr = collect( keys( ligdist ) )
   i = 1 # order iterator
   for ord in sortperm( collect( values( ligdist ) ) )
      outkey = keysarr[ord]
      print_heatrow( STDOUT, output[outkey]..., annobox[outkey], order="\t" * string(i) )
      i += 1
   end
end
###################################################################
main()
