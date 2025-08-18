# Use Ubuntu as base image
FROM ubuntu:22.04

# Avoid prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package list and install essential build tools
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    netcat-openbsd \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy source files (adjust as needed for your project structure)
COPY persistent_cell.cpp .
#COPY persistent_cell /app/persistent_cell
#COPY results.csv /app/results.csv

RUN g++ /app/persistent_cell.cpp -O3 -fomit-frame-pointer -fopenmp -m64 -std=c++11 -o /app/persistent_cell
#RUN chmod +x /app/persistent_cell
#RUN chmod +w /app/results.csv
RUN chmod -R 777 /app
#RUN touch /app/results.csv
#RUN chmod -R 777 /app/results.csv

EXPOSE 8080

# Compile the C++ program

# Make the executable file executable (if needed)
RUN chmod +x persistent_cell

# Run the program
#CMD ["/app/persistent_cell > /app/results.csv"]
#CMD ["/app/persistent_cell > results.csv"]
#CMD ["/app/persistent_cell"]
#CMD ["/bin/sh", "-c", "./persistent_cell > results.csv"]
#CMD ["/bin/sh", "-c", "./persistent_cell > results.csv && while true; do nc -l -p 8080 < results.csv; done"]
CMD ["/bin/sh", "-c", "./persistent_cell"]
#CMD ["sh", "-c", "./myprogram > output.txt 2>&1 && while true; do nc -l -p 8080 < output.txt; done"] 
