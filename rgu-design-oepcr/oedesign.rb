require 'bio'
require 'csv'
require "dnabitrue"
require 'fileutils'

include Dnabitrue

class MYNA < Bio::Sequence::NA
 # Calculate the melting temperature
 #  reference: http://www.sigmaaldrich.com/japan/lifescience/custom-products/custom-dna/oligo-learning-center.html#o04
 # 
 # s = Bio::Sequence::NA.new('CGACGTTGTAAAACGACGGCCAGT')
 # puts s.tm                            #=> 73.15
 # 
 # * (optional) _ct : conc of oligo nucleotide(mM) (default 0.5)
 # * (optional) _na : conc of Na+ (mM) (default 50)
 # * (optional) _fa : conc of Formamide (mol/L) (default 0)
 # *Returns*:: Float
 def tm(_ct=0.5, _na=50, _fa=0)
   naseq = self.dna
   len = naseq.length.to_f # base length
   ct = _ct.to_f/1000000   # concentration of oligo (mol/L)
   na = _na.to_f/1000      # concentration of Na+ (mol/L)
   fa = _fa.to_f           # concentration of Formamide (mol/L)
   @@r=1.987                   # gas constant R
   
   # to_upper
   naseq.upcase!

   # thermodynamic parameter 
   #  _h ...delta enthalpy
   #  _s ...delta entropy
   _h = {"AA" => -9.1, "TT" => -9.1, "AT" => -8.6, "TA" => -6.0,
            "CA" => -5.8, "TG" => -5.8, "GT" => -6.5, "AC" => -6.5,
            "CT" => -7.8, "AG" => -7.8, "GA" => -5.6, "TC" => -5.6,
            "CG" => -11.9, "GC" => -11.1, "GG" => -11.0, "CC" => -11.0}
   _s={"AA" => -24.0, "TT" => -24.0, "AT" => -23.9, "TA" => -16.9,
          "CA" => -12.9, "TG" => -12.9, "GT" => -17.3, "AC" => -17.3,
          "CT" => -20.8, "AG" => -20.8, "GA" => -13.5, "TC" => -13.5,
          "CG" => -27.8, "GC" => -26.7, "GG" => -26.6, "CC" => -26.6}
   tot_h=0.0
   tot_s=0.0
   for i in 0..(naseq.length-2)
     tot_h += _h[ naseq.slice(i,2) ]
     tot_s += _s[ naseq.slice(i,2) ]
   end

   count = self.composition
   at = count['a'] + count['t']
   gc = count['g'] + count['c']

   _tm = ((1000*tot_h)/(-10.8+tot_s+@@r*Math::log(ct/4)))-273.15+16.6*Math::log10(na)
   if  _tm > 80
     # gc% method
     tm = 81.5+16.6 * Math::log10(na) + 41*((gc.to_f)/len)-500/len-0.62*fa
   elsif _tm < 20
     # wallace method
     tm = 2*(at)+4*(gc)
   else
     tm = _tm
   end

   return tm
 end
end

def design_oeoligos(tseq)

	slen= tseq.length
	# (100..500).each{|slen|

		i=0
		ol=100
		while ol>65 do
			i=i+2
			ol=(slen+(19*(i-1)))/i
			# p ol	
		end
		out_ol=ol*i-slen
		add_len=out_ol/i
		mod_ol=ol

		p "#{slen}  #{i}  #{ol}  #{out_ol}   #{add_len}   #{mod_ol}"
		h_oligo=Hash.new

		tj=1
		h_oligo=Hash.new
		(1..i).each{|j|
			sj=mod_ol*(j-1)+1-18*(j-1)
			if j<i then
				ej=sj+mod_ol
				p "#{j}  #{sj}  #{ej}  #{ej-sj+1}" 

			else
				ej=slen
				p "#{j}  #{sj}  #{ej}  #{ej-sj+1}"
			end
			h_oligo[j]=[sj,ej]
			# p j_oligo
		}

	h_range2=Hash.new
	stj=0
	endj=0
	h_oligo.each{|j,jregion|
		# p jregion
		if j==1 then
			stoligo=jregion[0]
			endoligo=jregion[1]
			stj=endoligo-19
			endj=endoligo
		else
			stoligo=stj
			endoligo=h_range2[j][1]
			stj=endoligo-19
			endj=endoligo
		end
		# p "#{stoligo}, #{endoligo}, #{stj}, #{endj}"
		if j<i then
			jlen=19
			dtm=40
			lltm=52
			ultm=58  # some setting make infiinite roop
			
			while dtm < lltm || dtm > ultm do
				joligo= tseq[stj..endj]
				# p joligo
				# dna=Bio::Sequence::NA.new(joligo)
				# p dna.gc_percent
				dna = MYNA.new(joligo)
				dtm=dna.tm
				# p "#{joligo}  #{joligo.length}  #{dna.tm}"
				if dtm < lltm then
					jlen += 1
					if jlen.odd? then
						stj=stj-1
					else
						endj=endj+1
					end
				elsif dtm > ultm then
					jlen = jlen-1
					if jlen.odd? then
						stj=stj+1
					else
						endj=endj-1
					end
				end
			end
			# p "#{stj}, #{endj}"
			h_range2[j]=[stoligo,endj]
			h_range2[j+1]=[stj,h_oligo[j+1][1]]
		end
	}
	h_oligo2=Hash.new
	h_oligorange2=Hash.new
	h_range2.each{|j,olrange|
		# p "#{j}  #{olrange}  #{olrange[1]-olrange[0]+1}"
		h_oligo2[j]=tseq[olrange[0]-1..olrange[1]-1]
		h_oligorange2[j]=[olrange[0]-1,olrange[1]-1]
	}
	return h_oligo2,h_oligorange2

end


########################################################
# CSV.foreach("selected_primers_1.csv") do |row|
# 	p row
# end



amplfile=Dir.glob("tgtseq/*.fasta")[0]
p amplfile
basename=amplfile.sub("tgtseq/","").sub(".fasta","")
ff=Bio::FlatFile.auto(amplfile)
tseq=""


h_oeoligo=Hash.new
h_oeoligorange=Hash.new
h_oeoligo[basename]=Hash.new
h_oeoligorange[basename]=Hash.new
h_tgtseq=Hash.new
h_tgtseq[basename]=Hash.new
ff.each_entry do |f|
	sdef=f.definition
	sno=sdef
	# if sno==4 then
		
		p sno
		tseq=f.naseq
		p tseq
		
		h_oeoligo[basename][sno], h_oeoligorange[basename][sno]=design_oeoligos(tseq)

		h_tgtseq[basename][sno]=tseq.upcase
	# end
end
# p h_oeoligo

joinoligo=Hash.new
h_oeoligorange.each{|basename,data1|
	joinoligo[basename]=Hash.new
	data1.each{|sno,data2|
		
		jarray=Array.new
		data2.each{|olno,orange|
			jarray << orange.sort
		}
		jflat= jarray.flatten
		jflat.shift
		p jflat
		jflat.pop
		joinoligo[basename][sno]= jflat.each_slice(2).to_a.map{|x| x.sort}
	}
}

joligoseq=Hash.new
joinoligo.each{|basename,data1|
	joligoseq[basename]=Hash.new
	data1.each{|sno,jrange|
		joligoseq[basename][sno]=Array.new
		jrange.each_with_index{|rg,idx|
			joligoseq[basename][sno]<< [idx+1,h_tgtseq[basename][sno][rg[0]..rg[1]]]
		} 
	}
}

hitMode = 1
hitRate = 80
obj = FindPattern.new
obj.setHitMode( hitMode ) # hitMode　1: 曖昧　1以外：完全一致
obj.setHitRate( hitRate ) # hitrate  一致率

homol_join=Hash.new
joligoseq.each{|basename,data1|
	homol_join[basename]=Hash.new
	data1.each{|sno,data2|
		homol_join[basename][sno]=Array.new
		data2.each{|jo1|
			data2.each{|jo2|
				if jo1[0] < jo2[0] then
					flg=0
					# p "#{jo1[0]} #{jo2[0]}"
					j1len=jo1[1].length
					tgtlen=j1len/2
					p1jo=jo1[1][0,tgtlen]
					p2jo=jo1[1][tgtlen-1..-1]
					
					# p obj.FindPattern(p1jo,jo2[1])[0]
					# p obj.FindPattern(p2jo,jo2[1])[0]
					if obj.FindPattern(p1jo,jo2[1])[0] || obj.FindPattern(p2jo,jo2[1])[0] then
						flg=1
					end
					if flg==1 then 
						homol_join[basename][sno] << [jo1[0],jo2[0]]
					end
				end
			}
		}
	}
}

p homol_join
h_oeoligorange.each{|basename,data1|
	data1.each{|sno,data2|
		data2.each{|olno,rg|
			p "#{olno} #{rg[1]-rg[0]}"
		}
	}
}


h_oeol2=Hash.new
h_oeoligo.each{|tgv,data1|
	data1.each{|sno,data2|
		data2.each{|olno,oseq|
			odef="#{sno}|#{olno}"
			if olno.even? then
				onseq=Bio::Sequence::NA.new(oseq)
				oseq2=onseq.complement.to_s.upcase
				h_oeol2[odef]=oseq2
			else
				h_oeol2[odef]=oseq.upcase
			end
		}
	}
}

h_oeol2.each{|odef1,oseq1|
	h_oeol2.each{|odef2,oseq2|
		unless odef1==odef2 then
			if oseq1==oseq2 then
				p "#{odef1}  #{odef2}"
			end
		end
	}
}

CSV.open("order/order_oligo_#{basename}.csv","wb") do |csv|

	h_oeol2.each{|odef,oseq|
		csv << [odef,oseq,oseq.length]
	}

end

CSV.open("order/check_order_oligo_#{basename}.csv","wb") do |csv|
	presno=0

	h_oeol2.each{|odef,oseq|
		tgtsno=odef.split("|")[1].to_i
		unless tgtsno==presno then
			csv << ["Whole #{tgtsno}",h_tgtseq[basename][tgtsno]]
		end
		csv << [odef, oseq, oseq.length]
		presno=tgtsno
	}

end

homol_join.each{|basename,data1|
	CSV.open("order/homologous_oligo_#{basename}.csv","wb") do |csv|
	
		data1.each{|sno,data2|
			data2.each{|olno|
				csv << [basename,sno,olno.join(",")].flatten
			}
		}
	end
}

