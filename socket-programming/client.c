#include <winsock2.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#pragma comment(lib, "ws2_32.lib") // 自動鏈接 Winsock 庫

// 全域變數
volatile int Running = 1;
SOCKET clientSocket;

// 日誌函數
void log_message(const char *message) {
    printf("[LOG] %s\n", message);
}

// 錯誤處理函數
void handle_error(const char *errorMessage, int exitCode) {
    fprintf(stderr, "[ERROR] %s. Error Code: %d\n", errorMessage, WSAGetLastError());
    if (exitCode >= 0) {
        WSACleanup();
        exit(exitCode);
    }
}

// 初始化客戶端
SOCKET initialize_client(const char *serverIP, int port) {
    WSADATA wsa;
    SOCKET sock;
    struct sockaddr_in server;

    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) {
        handle_error("Failed to initialize Winsock", 1);
    }

    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == INVALID_SOCKET) {
        handle_error("Failed to create socket", 1);
    }

    server.sin_family = AF_INET;
    server.sin_addr.s_addr = inet_addr(serverIP);
    server.sin_port = htons(port);

    if (connect(sock, (struct sockaddr *)&server, sizeof(server)) == SOCKET_ERROR) {
        handle_error("Failed to connect to server", 1);
    }

    log_message("Successfully connected to server.");
    return sock;
}

// 處理用戶輸入的線程函數
DWORD WINAPI handle_input(LPVOID arg) {
    char userInput[256];

    while (Running) {
        fgets(userInput, sizeof(userInput), stdin);
        userInput[strcspn(userInput, "\n")] = '\0'; // 移除換行符

        if (strcmp(userInput, "end") == 0) {
            Running = 0;
            shutdown(clientSocket, SD_BOTH);
            closesocket(clientSocket);
            log_message("Closed the connection.");
            break;
        }

        if (send(clientSocket, userInput, strlen(userInput), 0) == SOCKET_ERROR) {
            handle_error("Failed to send message", -1);
        }
    }
    return 0;
}

// 處理伺服器消息
void handle_server_messages(SOCKET sock) {
    char buffer[256];
    int readSize;

    while (Running) {
        readSize = recv(sock, buffer, sizeof(buffer) - 1, 0);

        if (readSize > 0) {
            buffer[readSize] = '\0';
            printf("[SERVER]: %s", buffer);
        } else if (readSize == 0) {
            log_message("Server closed the connection.");
            Running = 0;
            break;
        }
    }
}

// 主程式
int main(void) {
    const char *serverIP = "127.0.0.1";
    int port = 5678;
    char clientName[256];
    HANDLE inputThread;

    clientSocket = initialize_client(serverIP, port);

    // 接收伺服器連線結果
    char serverResponse[256];
    int responseSize = recv(clientSocket, serverResponse, sizeof(serverResponse) - 1, 0);
    if (responseSize > 0) {
        serverResponse[responseSize] = '\0';
        if (strcmp(serverResponse, "full") == 0) {
            log_message("Maximum client limit reached. Connection refused.");
            closesocket(clientSocket);
            WSACleanup();
            return 1;
        }
    } else {
        handle_error("Error reading server response", 1);
    }

    // 輸入使用者名稱
    printf("Enter your name: ");
    fgets(clientName, sizeof(clientName), stdin);
    clientName[strcspn(clientName, "\n")] = '\0';

    if (send(clientSocket, clientName, strlen(clientName), 0) == SOCKET_ERROR) {
        handle_error("Failed to send name to server", 1);
    }

    // 啟動用戶輸入線程
    inputThread = CreateThread(NULL, 0, handle_input, NULL, 0, NULL);
    if (inputThread == NULL) {
        handle_error("Failed to create input thread", 1);
    }

    // 處理伺服器消息
    handle_server_messages(clientSocket);

    // 等待輸入線程結束
    WaitForSingleObject(inputThread, INFINITE);

    // 清理資源
    closesocket(clientSocket);
    WSACleanup();
    log_message("Client shut down.");
    return 0;
}
