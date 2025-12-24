CXX = g++
CXXFLAGS = -std=c++17 -I/usr/include/eigen3 -Wall -g -MMD -MP
TARGET = ns_fda_3D
SRCS = src/main.cpp src/inout.cpp src/mesh_n_model.cpp src/compute_flow.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include $(OBJS:.o=.d)

.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)