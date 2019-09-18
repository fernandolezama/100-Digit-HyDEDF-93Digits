function [ func ] = callFunction( funcNum )

funcs(1).name = 'Storns Chebyshev';        funcs(1).lowerlimit = -8192;     funcs(1).upperlimit = 8192;       funcs(1).dim = 9;         funcs(1).code=funcNum;
funcs(2).name = 'Inverse Hilbert';         funcs(2).lowerlimit = -16384;    funcs(2).upperlimit = 16384;      funcs(2).dim = 16;        funcs(2).code=funcNum;
funcs(3).name = 'Lennard-Jones';           funcs(3).lowerlimit = -4;        funcs(3).upperlimit = 4;          funcs(3).dim = 18;        funcs(3).code=funcNum;
funcs(4).name = 'Rastrigins';              funcs(4).lowerlimit = -100;      funcs(4).upperlimit = 100;        funcs(4).dim = 10;        funcs(4).code=funcNum;
funcs(5).name = 'Griewangks';              funcs(5).lowerlimit = -100;      funcs(5).upperlimit = 100;        funcs(5).dim = 10;        funcs(5).code=funcNum;
funcs(6).name = 'Weierstrass';             funcs(6).lowerlimit = -100;      funcs(6).upperlimit = 100;        funcs(6).dim = 10;        funcs(6).code=funcNum;
funcs(7).name = 'Modified Schwefels';      funcs(7).lowerlimit = -100;      funcs(7).upperlimit = 100;        funcs(7).dim = 10;        funcs(7).code=funcNum;
funcs(8).name = 'Expanded Schaffers F6';   funcs(8).lowerlimit = -100;      funcs(8).upperlimit = 100;        funcs(8).dim = 10;        funcs(8).code=funcNum;
funcs(9).name = 'Happy Cat';               funcs(9).lowerlimit = -100;      funcs(9).upperlimit = 100;        funcs(9).dim = 10;        funcs(9).code=funcNum;
funcs(10).name = 'Ackley';                 funcs(10).lowerlimit = -100;     funcs(10).upperlimit = 100;       funcs(10).dim = 10;       funcs(10).code=funcNum;

func = funcs(:,funcNum);

end
