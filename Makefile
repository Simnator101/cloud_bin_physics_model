IDIR = ./include
SDIR = ./src
BDIR = ./build

CC = gcc
CFLAGS = -g -Wall -std=c11
LIBS = -lm -lnetcdf

SRCS := $(wildcard $(SDIR)/*.c)
OBJS := $(patsubst $(SDIR)/%.c,$(BDIR)/%.o,$(SRCS))



all: SCModel

SCModel: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(BDIR)/%.o: $(SDIR)/%.c 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BDIR)/*.o 
	rm SCModel
