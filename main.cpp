#define _POSIX_SOURCE

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <semaphore.h>
#include <fcntl.h>

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <sstream>

char startData = '#', stopData = '*';
std::string cognatorPrefix = "COG";
std::string simulatorPrefix = "SIM";

struct IPCkit
{
   sem_t * semaphoreSim;
   sem_t * semaphoreCog;
   char * sharedMemory;
};

IPCkit init(key_t memoryKey, const char * semaphoreSimName, const char * semaphoreCogName, size_t memorySize)
{
   IPCkit ipcKit;

   //Getting shared memory
   int memoryId = shmget(memoryKey, memorySize, IPC_CREAT | 0666);
   if (memoryId < 0)
   {
      std::cout << "Error getting shared memory id" << std::endl;
      exit(1);
   }
   //Attached shared memory
   ipcKit.sharedMemory = (char *) shmat(memoryId, NULL, 0);
   if (ipcKit.sharedMemory == (char *) -1)
   {
      std::cout << "Error attaching shared memory id" << std::endl;
      exit(1);
   }

   ipcKit.semaphoreSim = sem_open(semaphoreSimName, O_CREAT, S_IRWXU, 1);
   ipcKit.semaphoreCog = sem_open(semaphoreCogName, O_CREAT, S_IRWXU, 0);

   return ipcKit;
}

void term(IPCkit ipcKit)
{
   sem_close(ipcKit.semaphoreCog);
   sem_close(ipcKit.semaphoreSim);
   shmdt(ipcKit.sharedMemory);
}

IPCkit get(key_t memoryKey, char * semaphoreSimName, char * semaphoreCogName, size_t memorySize)
{
   IPCkit ipcKit;

   //Getting shared memory
   int memoryId = shmget(memoryKey, memorySize, IPC_CREAT | 0666);
   if (memoryId < 0)
   {
      std::cout << "Error getting shared memory id" << std::endl;
      exit(1);
   }
   //Attached shared memory
   ipcKit.sharedMemory = (char *) shmat(memoryId, NULL, 0);
   if (ipcKit.sharedMemory == (char *) -1)
   {
      std::cout << "Error attaching shared memory id" << std::endl;
      exit(1);
   }

   ipcKit.semaphoreSim = sem_open(semaphoreSimName, O_CREAT);
   ipcKit.semaphoreCog = sem_open(semaphoreCogName, O_CREAT);

   return ipcKit;
}

void simulator(IPCkit ipcKit, unsigned long int maxInterval)
{
   //Loop
   unsigned long int X = 0;
   //Shared memory I/O ports
   char * sharedMemory = ipcKit.sharedMemory, sharedMemoryChar;
   std::stringstream readWriteOn;

   sem_wait(ipcKit.semaphoreSim);
   readWriteOn.str("");
   readWriteOn << startData << "Hello World" << X << stopData;
   memcpy(ipcKit.sharedMemory, readWriteOn.str().data(), readWriteOn.str().size());
   sem_post(ipcKit.semaphoreCog);

   while (X < maxInterval)
   {
      sem_wait(ipcKit.semaphoreSim);

      //Reading from the shared memory segment.
      sharedMemoryChar = ' ';
      readWriteOn.str("");
      sharedMemory = ipcKit.sharedMemory;
      do
      {
         memcpy(&sharedMemoryChar, sharedMemory++, 1);
         readWriteOn << sharedMemoryChar;
      }
      while (sharedMemoryChar != stopData);

      //Write on the shared memory segment
      std::cout << readWriteOn.str() << std::endl;
      readWriteOn.str("");
      readWriteOn << startData << "Hello World" << X << stopData;
      memcpy(ipcKit.sharedMemory, readWriteOn.str().data(), readWriteOn.str().size());

      sem_post(ipcKit.semaphoreCog);
      X++;
   }
}

void cognator(IPCkit ipcKit, unsigned long int maxInterval)
{
   //Loop
   unsigned long int X = 0;
   //Shared memory I/O ports
   char * sharedMemory = ipcKit.sharedMemory, sharedMemoryChar;
   std::stringstream readWriteOn;

   while (X < maxInterval)
   {
      sem_wait(ipcKit.semaphoreCog);

      //Reading from the shared memory segment.
      sharedMemoryChar = ' ';
      readWriteOn.str("");
      sharedMemory = ipcKit.sharedMemory;
      do
      {
         memcpy(&sharedMemoryChar, sharedMemory++, 1);
         readWriteOn << sharedMemoryChar;
      }
      while (sharedMemoryChar != stopData);

      //Write on the shared memory segment
      std::cout << readWriteOn.str() << std::endl;
      readWriteOn.str("");
      readWriteOn << startData << "Bye World" << X << stopData;
      memcpy(ipcKit.sharedMemory, readWriteOn.str().data(), readWriteOn.str().size());

      sem_post(ipcKit.semaphoreSim);
      X++;
   }
}

int main(int argc, const char **argv)
{
   key_t key;
   std::stringstream semaphoreSimName, semaphoreCogName;
   size_t memorySize;
   IPCkit ipcKit;
   unsigned long int X;
   unsigned short int sizeOfKey;

   if (argc == 5)
   {
      sizeOfKey = sizeof(argv[1])/sizeof(argv[1][0]);

      semaphoreSimName << simulatorPrefix;
      semaphoreCogName << cognatorPrefix;
      for (unsigned short int zeros = 0; zeros < sizeOfKey; zeros++)
      {
         semaphoreSimName << "0";
         semaphoreCogName << "0";
      }

      key = atoi(argv[1]);
      semaphoreSimName << argv[1];
      semaphoreCogName << argv[1];
      memorySize = atoi(argv[2]);
      X = atoi(argv[3]);
      ipcKit = init(key, semaphoreSimName.str().data(), semaphoreCogName.str().data(), memorySize);

      if (cognatorPrefix.compare(argv[4]) == 0) cognator(ipcKit, X);
      else if (simulatorPrefix.compare(argv[4]) == 0) simulator(ipcKit, X);

      term(ipcKit);
      sem_unlink(semaphoreSimName.str().data());
      sem_unlink(semaphoreCogName.str().data());
   }

   return 0;
}
