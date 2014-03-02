
# main compilation shell script:
g++ -I/Users/reynolds/develop/BAGEL/ ericompute.cc input_reader.cc /Users/reynolds/develop/BAGEL/src/integral/comprys/_co* -std=c++11 -g -o obj/compute

# also compiles bagel_interface.cc without also compiling Bagel - use for first round of debugging:
#g++ -I/Users/reynolds/develop/BAGEL/Debug -I/Users/reynolds/develop/BAGEL/ ericompute.cc input_reader.cc bagel_interface.cc /Users/reynolds/develop/BAGEL/src/integral/comprys/_co* -std=c++11 -g -o obj/compute
