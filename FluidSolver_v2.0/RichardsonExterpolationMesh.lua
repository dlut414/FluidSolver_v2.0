-- Richardson Exterpolation of Nu

print(" input: 40, 80, 160 \n");
print(" lambda: 2 \n");

print(" Choose the wall to exterpolate: \n");
print(" 0 : hot \n");
print(" 1 : cold \n");

str = io.read("*line");


path1 = "./40_upwind_Nu" .. str .. ".dat";
path2 = "./80_upwind_Nu" .. str .. ".dat";
path3 = "./160_upwind_Nu" .. str .. ".dat";
pathOut = "./RE_upwind_Nu" .. str .. ".dat";

h1 = (1./40); h2 = (1./80); h3 = (1./160); lambda = h1/h2;
X1_c, X2_c, X3_c = 0., 0., 0.; Xt_c = 0.;
X1_h, X2_h, X3_h = 0., 0., 0.; Xt_h = 0.;
K_c = 0.; K_h = 0.;
n_c = 0.; n_h = 0.;
C_c = 0.; C_h = 0.;


input1 = io.open(path1);
input2 = io.open(path2);
input3 = io.open(path3);
outputFile = io.open(pathOut, "w");

if(input1==nil or input2==nil or input3==nil) then
	print(" No file ", "\n");
else
	input1:read("*line");
	input2:read("*line");
	input3:read("*line");
	outputFile:write(" x y ux uy \n");
	while(true) do
		x1, y1, X1_h, X1_c = input1:read("*number", "*number", "*number", "*number");
		if(x1 == nil) then break end
		while(true) do
			x2, y2, X2_h, X2_c = input2:read("*number", "*number", "*number", "*number");
			if(math.abs(y1-y2) < 1e-5) then break end
		end
		while(true) do
			x3, y3, X3_h, X3_c = input3:read("*number", "*number", "*number", "*number");
			if(math.abs(y1-y3) < 1e-5) then break end
		end
		K_h = (X1_h - X2_h) / (X2_h - X3_h);
		K_c = (X1_c - X2_c) / (X2_c - X3_c);
		n_h = math.log(K_h) / math.log(lambda);
		n_c = math.log(K_c) / math.log(lambda);
		C_h = (X2_h - X3_h) / ((h3^n_h) - (h2^n_h));
		C_c = (X2_c - X3_c) / ((h3^n_c) - (h2^n_c));
		Xt_h = X3_h + C_h*(h3^n_h);
		Xt_c = X3_c + C_c*(h3^n_c);

		outputFile:write(string.format("%e %e %e %e\n", x1, y1, Xt_h, Xt_c));

		print(" n_h=", n_h, "\t", "n_c=", n_c, "\n");
		print(" C_h=", C_h, "\t", "C_c=", C_c, "\n");
		print(" Xt_h=", Xt_h, "\t", "Xt_c=", Xt_c, "\n");

	end
end

os.execute("pause");


