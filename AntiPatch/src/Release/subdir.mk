################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AntiPatch.cpp \
../Printer.cpp \
../Residue.cpp \
../SubPatcher.cpp \
../Table.cpp \
../main.cpp 

OBJS += \
./AntiPatch.o \
./Printer.o \
./Residue.o \
./SubPatcher.o \
./Table.o \
./main.o 

CPP_DEPS += \
./AntiPatch.d \
./Printer.d \
./Residue.d \
./SubPatcher.d \
./Table.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


