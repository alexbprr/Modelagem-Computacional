CC = g++
 
# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g -Wall -O3 
OTHER_FLAGS = -lm -lpthread -fopenmp -lboost_thread -lboost_system -lboost_timer -lboost_chrono -lboost_random 
# The build target 
TARGET = GillespieSolver
 
all: $(TARGET)

$(TARGET): $(TARGET).cpp
			$(CC) $(TARGET).cpp $(CFLAGS) -o $(TARGET) $(OTHER_FLAGS)
clean:
	$(RM) $(TARGET)