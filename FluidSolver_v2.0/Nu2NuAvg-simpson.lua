-- calculate Nu average

str = io.read("*line");

path0 = "./" .. str .. "_upwind_Nu0.dat";
path1 = "./" .. str .. "_upwind_Nu1.dat";

outputFile = io.open("./" .. str .. "_NuAvg.dat", "w");

outputFile:write("hot cold\n");

path = io.open(path0);
if(path==nil) then
	print(" No file 0\n")
else

	io.input(path)
		count = -1;
		m = 0;
		nu = 0;
		x = {};
		local pat = "(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s*"
		for n1, n2, n3, n4 in string.gfind(io.read("*all"), pat) do
			if(count ~= -1) then
				x[count] = n3;
			end
			count = count + 1;
		end
		for i=0,count-2,2 do
			nu = nu + x[i] + 4.0* x[i+1] + x[i+2];
		end
		nuAvg = nu/(3.0*(count-1));
		outputFile:write(nuAvg .. " ");
		print(nuAvg);
	io.close(path)

end

path = io.open(path1);
if(path==nil) then
	print(" No file 1\n")
else

	io.input(path)
		count = -1;
		m = 0;
		nu = 0;
		x = {};
		local pat = "(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s*"
		for n1, n2, n3, n4 in string.gfind(io.read("*all"), pat) do
			if(count ~= -1) then
				x[count] = n3;
			end
			count = count + 1;
		end
		for i=0,count-2,2 do
			nu = nu + x[i] + 4.0* x[i+1] + x[i+2];
		end
		nuAvg = nu/(3.0*(count-1));
		outputFile:write(nuAvg .. " ");
		print(nuAvg);
	io.close(path)

end
