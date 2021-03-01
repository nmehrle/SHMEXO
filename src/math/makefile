#
#  Copyright (c) 2020 by Cheng Li
#
AR := ar crl
CFLAGS := -O3
RANLIB := ranlib
LIB_FILE := libclimath.a
SRC_FILES := $(wildcard *.c)
OBJ_DIR := obj/
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(SRC_FILES:.c=.o)))

.PHONY : all clean $(OBJ_DIR)

all : $(OBJ_DIR) $(LIB_FILE)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB_FILE) : $(OBJ_FILES)
	$(AR) $(LIB_FILE) $(OBJ_FILES) ; \
	$(RANLIB) $(LIB_FILE)

clean :
	rm -rf $(OBJ_DIR)
	rm -f $(LIB_FILE)
