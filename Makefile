# executable name
TARGET = inference

# compiler
CXX = g++
CXXFLAGS = -std=c++14 -Wall

# source files
SRCS = source/inference.cpp source/functions.cpp source/objects.cpp

# objects files
OBJS = $(SRCS:.cpp=.o)

# required library
LIBS = -ltbb

# construction rules
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean
