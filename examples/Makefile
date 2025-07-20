# ========================================================================= #
#	                                                                    #
#	Makefile for building: Mathieu_ode_example                          #
#                                                                           #
# ========================================================================= #


SHELL = /bin/sh

VERSION = date "+%Y.%m%d%" 

DIR_C := ./src
DIR_N := ..
DIR_O := ./obj
 
CC      = gcc
DEBUG   = -Wall
CFLAGS  = -g -O0  
CFLAGS += $(DEBUG)   
LFLAGS  =   -l m  

TARGET  = Mathieu_ode_example  

$(DIR_O)/%.o : $(DIR_N)/%.c
	$(CC) $(CFLAGS) -c  $< -o   $@  

$(DIR_O)/%.o : $(DIR_C)/%.c
	$(CC) $(CFLAGS) -c  $< -o   $@  

$(TARGET) : $(DIR_O)/Mathieu_ode_example.o $(DIR_O)/HPGutil.o $(DIR_O)/NRutil.o  $(DIR_O)/HPGsignal.o $(DIR_O)/HPGmatrix.o   $(DIR_O)/HPGode.o
	$(CC) $(CFLAGS)  $^ -o   $@  $(LFLAGS)

install:
	chown root $(TARGET); chmod u+s $(TARGET); mv $(TARGET) /usr/local/bin/.

clean:
	rm $(DIR_O)/*.o 

