-- write Nu

path1 = "./out";
path2 = "./";

str = io.read("*line");
inputFileName = path1 .. "/" .. str;
inputFile = io.open(inputFileName, "r");

outputFileName1 = path2 .. "/" .. "120_upwind_Nu0.dat";
outputFileName2 = path2 .. "/" .. "120_upwind_Nu1.dat";
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
			if(count < 121) then
				outputFile1:write(n1 .. " " .. n2 .. " " .. n3 .. " " .. n4 .. "\n");
			else
				outputFile2:write(n1 .. " " .. n2 .. " " .. n3 .. " " .. n4 .. "\n");
			end
			count = count + 1;
		end
	io.close(inputFile)

end
