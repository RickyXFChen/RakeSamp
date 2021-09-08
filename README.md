# RakeSamp
RakeSamp Software for sampling RNA secondary structures

To be specific, the program produces a uniformly random sample of RNA secondary structure with a prescribed number of hairpin loops, bulges, interior loops of distinct length, and a prescribed number of multiloops of distinct degree and length, and a prescribed number of helices of distinct size, and an exterior loop of a prescibed degree and length.

The only input is the file name of the .txt file where the corresponding parameters are specified. The format of these data in the .txt file is as follows:
1. The first row specifies the degree and the length of the exterior loop, separated by a space;
2. The second row specifies the number of helices of size i>=1, separated by a space;
3. The third row specifies the number of hairpins of length i>=0, separated by a space;
4. The fouth row specifies the number of bulges of length i>=0, separated by a space;
5. The fifth row specifies the number of interior loops of lenght i>=0, separated by a space;
6. Starting from the 6th row, the first entry of each row specifies the degree of the multiloops, then followed entries specify the number of multiloops of length i>=0 of that degree, separated by spaces.

Please refer to an example in the file exam.txt
