-- calculate Nu max

str = io.read("*line");

path0 = "./" .. str .. "_upwind_Nu0.dat";
path1 = "./" .. str .. "_upwind_Nu1.dat";

outputFile = io.open("./" .. str .. "_NuMax.dat", "w");

outputFile:write("hot cold\n");

path = io.open(path0);
if(path==nil) then
	print(" No file 0\n")
else

	io.input(path)
		count = -1;
		nu = 0;
		local pat = "(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s*"
		for n1, n2, n3, n4 in string.gfind(io.read("*all"), pat) do
			if(count ~= -1) then
				if(math.abs(tonumber(n3)) > math.abs(nu)) then
					nu = n3
				end
			end
			count = count + 1;
		end
		outputFile:write(nu .. " ");
	io.close(path)

end

path = io.open(path1);
if(path==nil) then
	print(" No file 1\n")
else

	io.input(path)
		count = -1;
		nu = 0;
		local pat = "(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s*"
		for n1, n2, n3, n4 in string.gfind(io.read("*all"), pat) do
			if(count ~= -1) then
				if(math.abs(tonumber(n3)) > math.abs(nu)) then
					nu = n3
				end
			end
			count = count + 1;
		end
		outputFile:write(nu .. " ");
	io.close(path)

end
