obj = main.cpp lib.cpp lib.h
cxx = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
gpr = -pg $(cxx)
wo3 = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
testfold = ../matrixtests/

a.out: $(obj)

	g++ $(cxx) $(obj) -o a.out

zip:

	zip Frolov_PS.zip $(obj) makefile

gp.out: $(obj)

	g++ $(gpr) $(obj) -o gp.out


ga.out: $(obj)

	g++ $(wo3) $(obj) -o ga.out

git:
	./commit.sh
	

push:
	git push

test:
	./test.sh  $(testfold)

