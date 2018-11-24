require 'csv'
require 'oj'


b_no=Hash.new
CSV.foreach("base_number.csv") do |row|
	b_no[row[0]]=row[1].to_i
end
p b_no
ncomp=Hash.new
CSV.foreach("complement_base.csv") do |row|
	ncomp[row[0].to_i]=row[1].to_i
end

cpp=Hash.new
CSV.foreach("complement_pair_point.csv") do |row|
	cpp[[row[0],row[1]]]=row[2].to_i
end
ncpp=Hash.new
cpp.each{|set,pt|
	ncpp[[b_no[set[0]],b_no[set[1]]]]=pt
}

templ=Hash.new
CSV.foreach("再合成Template.csv") do |row|
	templ[row[0]]=Hash.new
	templ[row[0]]["s"]=row[4].split("").map{|x| b_no[x]}
	templ[row[0]]["r"]=templ[row[0]]["s"].map{|x| ncomp[x]}
end

dspoint=Hash.new
ndspoint=Hash.new
templ.each{|tn,d1|
	dspoint[tn]=Hash.new
	ndspoint[tn]=Hash.new
	(0..d1["r"].length-10).to_a.each{|i|
		dspoint[tn][i]=Hash.new
		ps=d1["s"][i,10]
		(0..d1["r"].length-10).to_a.each{|j|
			ppoint=ps.zip(d1["r"][j,10]).map{|set| ncpp[set]}.inject(:+)
			if ppoint > 7 then
				dspoint[tn][i][j]=ppoint
			end
		}
		ndspoint[tn][i]=dspoint[tn][i].length
	}
}
ndspoint.each{|tn,d|
	p tn
	p d.select{|k,v| v==1}
}