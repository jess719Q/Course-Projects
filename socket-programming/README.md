# Simple TCP Chat in C (Windows)

This project implements a basic TCP chat system for Windows using **Winsock2**. It consists of a **server** and a **client** program. Multiple clients can connect to the server and exchange messages in real-time.

## Features

- Server can handle up to 10 clients simultaneously.
- Clients send their names upon connecting.
- Messages are broadcasted to all connected clients.
- Non-blocking sockets ensure smooth communication.
- Server and client both support shutting down gracefully by typing `end`.

## Usage

### Server
1. Compile the server program.
2. Run the executable.
3. Type `end` to shut down the server.

### Client
1. Compile the client program.
2. Run the executable.
3. Enter your name when prompted.
4. Type messages to send to all connected clients.
5. Type `end` to disconnect.
