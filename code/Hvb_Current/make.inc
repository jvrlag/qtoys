
### MIN. CONFIGURATIONS
# LIBRARY NAME
NAME=Hvb
# LIBRARY VERSION
VERSION=Current
# SOURCES PATH
SRCD=./
# OUTPUT PATH BASE
PTHBASE=/usr/local/
# COMPILER
CC=g++

### OPT. CONFIGURATIONS
## DEV FLAGS
#FLAGS=-pedantic -Wall -g -O0
## PROD FLAGS
FLAGS=-pedantic -Wall -O3

### EXPERT CONFS.
GNAME=lib$(NAME).a
LNAME=lib$(NAME)-v$(VERSION).a
LIBD=$(PTHBASE)lib/
INCD=$(PTHBASE)include/

HVBSRC=$(wildcard *.cc)
HVBBSE=$(basename $(HVBSRC))
HVBS=$(addsuffix .o,$(HVBBSE))
HVBH=$(addsuffix .h,$(HVBBSE))

