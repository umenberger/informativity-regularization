function M = ss2mat(sys)

    M = [sys.A, sys.B;
         sys.C, sys.D];

end