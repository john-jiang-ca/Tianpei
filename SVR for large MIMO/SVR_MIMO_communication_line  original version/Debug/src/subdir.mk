################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/SVR_MIMO_communication_line.cpp \
../src/svm.cpp 

C_SRCS += \
../src/svm-predict.c \
../src/svm-scale.c \
../src/svm-train.c 

OBJS += \
./src/SVR_MIMO_communication_line.o \
./src/svm-predict.o \
./src/svm-scale.o \
./src/svm-train.o \
./src/svm.o 

C_DEPS += \
./src/svm-predict.d \
./src/svm-scale.d \
./src/svm-train.d 

CPP_DEPS += \
./src/SVR_MIMO_communication_line.d \
./src/svm.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


