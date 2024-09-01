# Use an official Python runtime as a parent image
FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Install dependencies
RUN apt-get update && \
    apt-get install -y python3-tk && \
    apt-get install -y snapd && \
    snap install root-framework

# Install PyROOT
#RUN pip install pyroot

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Run the application
CMD ["pyroot", "./disp_GUI.py"]
