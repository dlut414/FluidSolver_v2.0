-- Richardson Exterpolation of Nu

print(" input: 80, 120, 180 \n");
print(" lambda: 1.5 \n");

print(" Choose the value to exterpolate: \n");
print(" 1 : NuAvg \n");
print(" 2 : NuMax \n");

item = io.read("*line");

str = nil;
if(item == "1") then
	str = "NuAvg";
elseif(item == "2") then
	str = "NuMax";
end


path1 = "./" .. "80_" .. str .. ".dat";
path2 = "./" .. "120_" .. str .. ".dat";
path3 = "./" .. "180_" .. str .. ".dat";

h1 = (1./80); h2 = (1./120); h3 = (1./180); lambda = h1/h2;
X1_c, X2_c, X3_c = 0., 0., 0.; Xt_c = 0.;
X1_h, X2_h, X3_h = 0., 0., 0.; Xt_h = 0.;
K_c = 0.; K_h = 0.;
n_c = 0.; n_h = 0.;
C_c = 0.; C_h = 0.;


path = io.open(path1);
if(path==nil) then
	print(" No file ", path1, "\n");
else
	io.input(path)
		io.read("*line");
		X1_h, X1_c = io.read("*number", "*number");
	io.close(path)
end
path = io.open(path2);
if(path==nil) then
	print(" No file ", path2, "\n");
else
	io.input(path)
		io.read("*line");
		X2_h, X2_c = io.read("*number", "*number");
	io.close(path)
end
path = io.open(path3);
if(path==nil) then
	print(" No file ", path3, "\n");
else
	io.input(path)
		io.read("*line");
		X3_h, X3_c = io.read("*number", "*number");
	io.close(path)
end

K_c = (X1_c - X2_c) / (X2_c - X3_c);
K_h = (X1_h - X2_h) / (X2_h - X3_h);
n_c = math.log(K_c) / math.log(lambda);
n_h = math.log(K_h) / math.log(lambda);
C_c = (X2_c - X3_c) / ((h3^n_c) - (h2^n_c));
C_h = (X2_h - X3_h) / ((h3^n_h) - (h2^n_h));
Xt_c = X3_c + C_c*(h3^n_c);
Xt_h = X3_h + C_h*(h3^n_h);

outputFile = io.open("./" .. str .. "_RE.dat", "w");
outputFile:write(" hot \t cold\n");
outputFile:write(" n_h=", n_h, "\t", "n_c=", n_c, "\n");
outputFile:write(" C_h=", C_h, "\t", "C_c=", C_c, "\n");
outputFile:write(" Xt_h=", Xt_h, "\t", "Xt_c=", Xt_c, "\n");



