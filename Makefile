# Makefile for profile

NAME	= profile
VERSION	= $(shell git describe --tags --long)

CC		= gcc
CFLAGS	= -O3 -mcmodel=medium -Wall -pedantic -fopenmp -I$(LOCAL_LIB_PATH)/include -DVERSION=\"${VERSION}\"
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -lgsl -lgslcblas -liof -lartsfc

SRCS	= $(wildcard *.c)

# Rules

$(NAME): $(SRCS) Makefile
	$(CC) $(CFLAGS) $(SRCS) -o $(NAME) $(LIBS)

clean:
	rm -f *~ $(NAME)
