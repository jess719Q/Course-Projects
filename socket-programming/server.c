#include <winsock2.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#pragma comment(lib, "ws2_32.lib") // 自動鏈接 Winsock 庫

#define MAX_CLIENTS 10
#define BUFFER_SIZE 256

volatile int serverRunning = 1; // 控制伺服器運行狀態
SOCKET clientSockets[MAX_CLIENTS];
BOOL clientActive[MAX_CLIENTS] = {0};

// 初始化 Winsock 和伺服器 Socket
SOCKET initialize_server(const char *ip, int port) {
    WSADATA wsa;
    SOCKET serverSocket;
    struct sockaddr_in serverAddr;

    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) {
        printf("WSAStartup failed. Error Code: %d\n", WSAGetLastError());
        exit(EXIT_FAILURE);
    }

    serverSocket = socket(AF_INET, SOCK_STREAM, 0);
    if (serverSocket == INVALID_SOCKET) {
        printf("Socket creation failed. Error Code: %d\n", WSAGetLastError());
        WSACleanup();
        exit(EXIT_FAILURE);
    }

    serverAddr.sin_family = AF_INET;
    serverAddr.sin_addr.s_addr = inet_addr(ip);
    serverAddr.sin_port = htons(port);

    if (bind(serverSocket, (struct sockaddr*)&serverAddr, sizeof(serverAddr)) == SOCKET_ERROR) {
        printf("Bind failed. Error Code: %d\n", WSAGetLastError());
        closesocket(serverSocket);
        WSACleanup();
        exit(EXIT_FAILURE);
    }

    // 設置非阻塞模式
    u_long mode = 1; // 1 表示非阻塞模式
    ioctlsocket(serverSocket, FIONBIO, &mode);

    if (listen(serverSocket, MAX_CLIENTS) == SOCKET_ERROR) {
        printf("Listen failed. Error Code: %d\n", WSAGetLastError());
        closesocket(serverSocket);
        WSACleanup();
        exit(EXIT_FAILURE);
    }

    printf("Server listening on %s:%d...\n", ip, port);
    return serverSocket;
}

// 廣播訊息給所有活躍客戶端
void broadcast_message(const char *message, SOCKET excludeSocket) {
    for (int i = 0; i < MAX_CLIENTS; i++) {
        if (clientActive[i] && clientSockets[i] != excludeSocket) {
            send(clientSockets[i], message, strlen(message), 0);
        }
    }
}

// 處理單個客戶端的通訊
DWORD WINAPI handle_client(LPVOID param) {
    int clientIndex = *(int*)param;
    SOCKET clientSocket = clientSockets[clientIndex];
    char buffer[BUFFER_SIZE];
    char clientName[BUFFER_SIZE];
    char message[BUFFER_SIZE];
    int bytesRead;

    // 接收客戶端名稱
    while(TRUE){
        bytesRead = recv(clientSocket, clientName, sizeof(clientName) - 1, 0);
        if (bytesRead == SOCKET_ERROR) {
                int error = WSAGetLastError();
                if (error == WSAEWOULDBLOCK) {
                    Sleep(100);
                    continue;
                } else {
                    printf("Recv error from client %d: %d\n", clientIndex, error);
                    break;
                }
        }else if (bytesRead == 0) {
            printf("Client %d disconnected during name exchange.\n", clientIndex);
            clientActive[clientIndex] = FALSE;
            closesocket(clientSocket);
            break;
        }
        
        break;
    }
    clientName[bytesRead] = '\0';

    snprintf(message, sizeof(message), "%s has joined the chat!\n", clientName);
    printf("%s", message);
    send(clientSocket, message, strlen(message), 0);
    broadcast_message(message, clientSocket);

    // 處理訊息傳遞
    while (serverRunning) {
        bytesRead = recv(clientSocket, buffer, sizeof(buffer) - 1, 0);
        if (bytesRead == SOCKET_ERROR) {
            int error = WSAGetLastError();
            if (error == WSAEWOULDBLOCK) {
                Sleep(100);
                continue;
            } else {
                printf("Recv error from client %d: %d\n", clientIndex, error);
                break;
            }
        } else if (bytesRead == 0) {
            // 客戶端斷開連接
            snprintf(message, sizeof(message), "%s has left the chat!\n", clientName);
            printf("%s", message);
            broadcast_message(message, clientSocket);
            break;
        }

        buffer[bytesRead] = '\0';
        snprintf(message, sizeof(message), "%s: %s\n", clientName, buffer);
        printf("%s", message);
        broadcast_message(message, clientSocket);
    }

    clientActive[clientIndex] = FALSE;
    closesocket(clientSocket);
    return 0;
}

// 處理伺服器輸入
DWORD WINAPI handle_input(LPVOID arg) {
    char input[BUFFER_SIZE];
    while (serverRunning) {
        fgets(input, sizeof(input), stdin);
        input[strcspn(input, "\n")] = '\0'; // 移除換行符
        if (strcmp(input, "end") == 0) {
            serverRunning = 0;
            printf("Server shutting down...\n");
            break;
        }
    }
    return 0;
}

int main() {
    SOCKET serverSocket = initialize_server("127.0.0.1", 5678);
    struct sockaddr_in clientAddr;
    int addrLen = sizeof(clientAddr);
    HANDLE clientThreads[MAX_CLIENTS];
    HANDLE inputThread;
    int clientIndex;

    // 啟動輸入線程
    inputThread = CreateThread(NULL, 0, handle_input, NULL, 0, NULL);

    // 接受客戶端連接
    while (serverRunning) {
        SOCKET clientSocket = accept(serverSocket, (struct sockaddr*)&clientAddr, &addrLen);
        if (clientSocket == INVALID_SOCKET) {
            if (WSAGetLastError() == WSAEWOULDBLOCK) {
                Sleep(100);
                continue;
            } else {
                printf("Accept failed. Error Code: %d\n", WSAGetLastError());
                break;
            }
        }

        // 分配客戶端到空閒插槽
        clientIndex = -1;
        for (int i = 0; i < MAX_CLIENTS; i++) {
            if (!clientActive[i]) {
                clientIndex = i;
                clientActive[i] = TRUE;
                clientSockets[i] = clientSocket;
                break;
            }
        }

        if (clientIndex == -1) {
            printf("Maximum client limit reached. Connection refused.\n");
            send(clientSocket, "full", strlen("full"), 0);
            closesocket(clientSocket);
        } else {
            printf("Client %d connected.\n", clientIndex);
            send(clientSocket, "avail", strlen("avail"), 0);
            CreateThread(NULL, 0, handle_client, &clientIndex, 0, NULL);
        }
    }

    WaitForSingleObject(inputThread, INFINITE);

    // 關閉所有客戶端
    for (int i = 0; i < MAX_CLIENTS; i++) {
        if (clientActive[i]) {
            closesocket(clientSockets[i]);
        }
    }

    closesocket(serverSocket);
    WSACleanup();
    printf("Server shut down.\n");
    return 0;
}
