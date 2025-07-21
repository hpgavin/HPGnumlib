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

all : Mathieu_ode_example  Hysteresis_ode_example  

$(DIR_O)/%.o : $(DIR_N)/%.c
	$(CC) $(CFLAGS) -c  $< -o   $@  

$(DIR_O)/%.o : $(DIR_C)/%.c
	$(CC) $(CFLAGS) -c  $< -o   $@  

Mathieu_ode_example : $(DIR_O)/Mathieu_ode_example.o $(DIR_O)/HPGutil.o $(DIR_O)/NRutil.o  $(DIR_O)/HPGsignal.o $(DIR_O)/HPGmatrix.o   $(DIR_O)/HPGode.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

Hysteresis_ode_example : $(DIR_O)/Hysteresis_ode_example.o  $(DIR_O)/HPGutil.o $(DIR_O)/NRutil.o  $(DIR_O)/HPGsignal.o $(DIR_O)/HPGmatrix.o   $(DIR_O)/HPGode.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

#	$(CC) $(CFLAGS) -c  $< -o $@  $(LFLAGS)
#	$(CC) $(CFLAGS) -c  $^ -o $@  $(LFLAGS)


install:
	chown root Mathieu_ode_example  Hysteresis_ode_example ; chmod u+s Mathieu_ode_example  Hysteresis_ode_example  ; mv Mathieu_ode_example  Hysteresis_ode_example   /usr/local/bin/.

clean:
	rm $(DIR_O)/*.o 

