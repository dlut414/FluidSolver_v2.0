-- write velocity

path1 = "./out";
path2 = "C:/Users/HU FANGYUAN/Documents/publish/texGraph/data/thermalCavity/Ra1E6/p2r0q2s2tIE_static_upwind";

str = io.read("*line");
inputFileName = path1 .. "/" .. str;
inputFile = io.open(inputFileName, "r");

outputFileName1 = path2 .. "/" .. "120_upwind_ux.dat";
outputFileName2 = path2 .. "/" .. "120_upwind_uy.dat";
outputFile1 = io.open(outputFileName1, "w");
outputFile2 = io.open(outputFileName2, "w");

outputFile1:write("x y ux uy\n");
outputFile2:write("x y ux uy\n");

if(inputFile==nil) then
	print(" No file \n")
else

	io.input(inputFile)
		count = 0;
		local pat = "(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s*"
		for n1, n2, n3, n4 in string.gfind(io.read("*all"), pat) do
			if(count < 101) then
				outputFile1:write(n1 .. " " .. n2 .. " " .. n3 .. " " .. n4 .. "\n");
			else
				outputFile2:write(n1 .. " " .. n2 .. " " .. n3 .. " " .. n4 .. "\n");
			end
			count = count + 1;
		end
	io.close(inputFile)

end
