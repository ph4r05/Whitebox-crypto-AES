CXXFLAGS =	-O0 -g -Wall -fmessage-length=0 -L -lntl -std=c++0x

OBJS =		MGR_NTL.o GenericAES.o NTLUtils.o MixingBijections.o WBAES.o WBAESGenerator.o

LIBS = -lntl

TARGET =	MGR_NTL

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
