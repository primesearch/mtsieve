/* main.c -- (C) Mark Rodenkirch, November 2017

   Single-threaded CPU sieve application framework.

   For each prime p in 3 <= p0 <= p < p1 < 2^62
     Do something with p
     
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <stdarg.h>
#ifndef WIN32
#include <sys/resource.h>
#endif
#include <time.h>
#include <signal.h>
#include "main.h"

#if defined (WIN32)
#include <Windows.h>
#if defined (_MSC_VER) && defined (MEMLEAK)
   _CrtMemState mem_dbg1;
   HANDLE hLogFile;
#endif
#endif

#include "App.h"

int   ProcessArgs(App *theApp, int argc, char *argv[]);
void  FatalError(const char *fmt, ...);
void  MemoryLeakEnter(void);
void  MemoryLeakExit(void);

// Local variables
bool  help_opt = false;

int64_t cpuBytes = 0;
App  *theApp = 0;

volatile bool gb_ForceQuit = false;

void SetQuitting(int sig)
{
   if (gb_ForceQuit)
   {
      fprintf(stderr, "\nCTRL-C already hit, so I'm going to terminate now.  Remaining terms file is not complete\n");
      exit(0);
   }
   
   theApp->Interrupt();
   gb_ForceQuit = true;
}

int   main(int argc, char *argv[])
{
   MemoryLeakEnter();

   cpuBytes = 0;

   theApp = get_app();

   theApp->Banner();
   
   ProcessArgs(theApp, argc, argv);

   // Ignore SIGHUP, as to not die on logout.
   // We log to file in most cases anyway.
#ifdef SIGHUP
   signal(SIGHUP, SIG_IGN);
#endif
#ifdef SIGHUP
   signal(SIGQUIT, SetQuitting);
#endif
   signal(SIGINT, SetQuitting);
   signal(SIGTERM, SetQuitting);

#ifdef _WIN32
   SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);
   SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_IDLE);
#else
   // This seems backwards, but PRIO_MAX refers to maximum niceness.
   setpriority(PRIO_PROCESS, 0, PRIO_MAX);
#endif

   theApp->Run();

   delete theApp;

   MemoryLeakExit();
   return 0;
}

static void help(void)
{
   printf("-h --help             prints this help\n");
}

static int parse_option(App *theApp, int opt, char *arg, const char *source)
{
  int status = 0;

  switch (opt)
  {
    case 'h':
      help_opt = true;
      break;

    default:
      status = theApp->ParseOption(opt, arg, source);
      break;
  }

  return status;
}

// Process command-line options using getopt_long().
// Non-option arguments are treated as if they belong to option zero.
// Returns the number of options processed.
int   ProcessArgs(App *theApp, int argc, char *argv[])
{
   int      count = 0, ind = -1, opt;
   string   shortOpts = "h";
   struct option longOpts[50];

   memset(longOpts, 0, 50*sizeof(struct option));

   AppendLongOpt(longOpts, "help",        no_argument,       0, 'h');

   theApp->AddCommandLineOptions(shortOpts, longOpts);

   while ((opt = getopt_long(argc, argv, shortOpts.c_str(), longOpts, &ind)) != -1)
      switch (parse_option(theApp, opt, optarg, NULL))
      {
         case 0:
            ind = -1;
            count++;
            break;

         case P_FAILURE:
            // If ind is unchanged then this is a short option, otherwise long.
            if (ind == -1)
               FatalError("%s: invalid argument -%c %s", argv[0], opt, optarg);
            else
               FatalError("%s: invalid argument --%s %s", argv[0], longOpts[ind].name, optarg);

         case P_OUT_OF_RANGE:
            // If ind is unchanged then this is a short option, otherwise long.
            if (ind == -1)
               FatalError("%s: out of range argument -%c %s", argv[0], opt, optarg);
            else
               FatalError("%s: out of range argument --%s %s", argv[0], longOpts[ind].name, optarg);

         default:
            FatalError("%s: invalid option %s", argv[0], argv[optind]);
      }

   while (optind < argc)
      switch (parse_option(theApp, 0, argv[optind], NULL))
      {
         case 0:
            optind++;
            count++;
            break;

         case P_FAILURE:
            FatalError("%s: invalid non-option argument %s", argv[0], argv[optind]);

         case P_OUT_OF_RANGE:
            FatalError("%s: out of range non-option argument %s", argv[0], argv[optind]);

         default:
            FatalError("%s: invalid option %s", argv[0], argv[optind]);
      }

   if (help_opt)
   {
      help();
      theApp->Help();
   }

   return count;
}

void     FatalError(const char *fmt, ...)
{
   va_list  args;
   char     buffer[20000];

   va_start(args, fmt);
   vsprintf(buffer, fmt, args);
   va_end(args);

   fprintf(stderr, "Fatal Error:  %s\n", buffer);
   exit(-1);
}

void  AppendLongOpt(struct option *longOpts, const char *name, int has_arg, int *flag, char charSwitch)
{
   while (longOpts->name != 0)
      longOpts++;

   longOpts->name = name;
   longOpts->has_arg = has_arg;
   longOpts->flag = flag;
   longOpts->val = charSwitch;
}

void *xmalloc(size_t requestedSize)
{
   void     *allocatedPtr;
   char     *currentPtr;
   size_t    allocatedSize = requestedSize + 150;
   uint64_t  temp;

   // Allocate extra memory because we need to align to a 64-byte boundary
   // as the memory might be referenced by AVX instructions.  We will also
   // put the original pointer for the allocated memory and the orignal size
   // for the allocated memory at the front of this and follow by 0xff to
   // verify that someone isn't running past the end of their allocated memory.
   if ((allocatedPtr = malloc(allocatedSize)) == NULL)
      FatalError("Unable to allocate %" PRId64" bytes of memory", (int64_t) requestedSize);

   cpuBytes += allocatedSize;

   memset(allocatedPtr, 0x00, allocatedSize);
   
   temp = (uint64_t) allocatedPtr;
   
   // We want temp divisible by 64 and within the allocated area
   temp = temp - (temp%64) + 64;
   
   // Align to a 64-byte boundary
   currentPtr = (char *) temp;

   // Put the pointer to the allocated memory here
   *(uint64_t *) currentPtr = (uint64_t) allocatedPtr;
      
   // Put the size of the  to what was actually allocated here
   *(size_t *) (currentPtr + 8) = allocatedSize;

   // Get to the next boundary
   currentPtr += 64;
   
   // This will help us detect someone going past what they are supposed to
   *(currentPtr + requestedSize) = 0xff;

   return (void *) currentPtr;
}

void xfree(void *memoryPtr)
{
   char     *currentPtr;
   void     *allocatedPtr;
   size_t    allocatedSize;
   uint64_t  temp;
   
   temp = (uint64_t) memoryPtr;
   
   currentPtr = (char *) (temp - 64);
   
   // Reduce by what we actually allocated
   allocatedSize = *(size_t *) (currentPtr + 8);
   cpuBytes -= allocatedSize;
   
   allocatedPtr =  (void *) *(uint64_t *) currentPtr;
   
   free(allocatedPtr);
}

void  MemoryLeakEnter(void)
{
   #if defined (WIN32) && defined (_MSC_VER) && defined (MEMLEAK)
   hLogFile = CreateFile("KBNCAppcl_memleak.txt", GENERIC_WRITE,
      FILE_SHARE_WRITE, NULL, CREATE_ALWAYS,
      FILE_ATTRIBUTE_NORMAL, NULL);

   _CrtSetReportMode( _CRT_WARN, _CRTDBG_MODE_FILE);
   _CrtSetReportFile( _CRT_WARN, hLogFile);
   _CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ERROR, hLogFile );
   _CrtSetReportMode( _CRT_ASSERT, _CRTDBG_MODE_FILE );
   _CrtSetReportFile( _CRT_ASSERT, hLogFile );

   // Store a memory checkpoint in the s1 memory-state structure
   _CrtMemCheckpoint( &mem_dbg1 );

   // Get the current state of the flag
   int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG );
   // Turn Off (AND) - DO NOT Enable debug heap allocations and use of memory block type identifiers such as _CLIENT_BLOCK
   //tmpFlag &= ~_CRTDBG_ALLOC_MEM_DF;
   tmpFlag |= _CRTDBG_ALLOC_MEM_DF;

   // Turn Off (AND) - prevent _CrtCheckMemory from being called at every allocation request.  _CurCheckMemory must be called explicitly
   tmpFlag &= ~_CRTDBG_CHECK_ALWAYS_DF;
   //tmpFlag |= _CRTDBG_CHECK_ALWAYS_DF;

   // Turn Off (AND) - Do NOT include __CRT_BLOCK types in leak detection and memory state differences
   tmpFlag &= ~_CRTDBG_CHECK_CRT_DF;
   //tmpFlag |= _CRTDBG_CHECK_CRT_DF;

   // Turn Off (AND) - DO NOT Keep freed memory blocks in the heap’s linked list and mark them as freed
   tmpFlag &= ~_CRTDBG_DELAY_FREE_MEM_DF;

   // Turn Off (AND) - Do NOT perform leak check at end of program run.
   //tmpFlag &= ~_CRTDBG_LEAK_CHECK_DF;
   tmpFlag |= _CRTDBG_LEAK_CHECK_DF;

   // Set the new state for the flag
   _CrtSetDbgFlag( tmpFlag );
#endif
}
void  MemoryLeakExit(void)
{
#if defined (WIN32) && defined (_MSC_VER) && defined (MEMLEAK)
   bool  deleteIt;
   _CrtMemState mem_dbg2, mem_dbg3;

   if ( _CrtMemDifference( &mem_dbg3, &mem_dbg1, &mem_dbg2 ) )
   {
      _RPT0(_CRT_WARN, "\nDump the changes that occurred between two memory checkpoints\n");
      _CrtMemDumpStatistics( &mem_dbg3 );
      _CrtDumpMemoryLeaks( );
      deleteIt = false;
   }
   else
      deleteIt = true;

   CloseHandle(hLogFile);

   if (deleteIt)
      unlink("prpclient_memleak.txt");
#endif
}
