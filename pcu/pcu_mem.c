/****************************************************************************** 

  (c) 2023 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include <PCU.h>

#if defined(__APPLE__)

#include <mach/task.h>
#include <mach/mach_init.h>

#elif defined(__bgq__)

//the BG/Q headers have warning-triggering
//code in them.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <spi/include/kernel/memory.h>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#else

#include <malloc.h> //warning - this is GNU-specific

#endif

double PCU_GetMem() {
  const double M = 1024*1024;
#if defined(__APPLE__)
  bool resident = true;
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
  size_t size = (resident ? t_info.resident_size : t_info.virtual_size);
  return (double)size/M;
#elif defined(__bgq__)
  size_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  return (double)heap/M;
#elif defined(__GNUC__) && defined(PUMI_HAS_MALLINFO2)
  struct mallinfo2 meminfo_now = mallinfo2();
  return ((double)meminfo_now.arena)/M;
#elif defined(__GNUC__)
  struct mallinfo meminfo_now = mallinfo();
  return ((double)meminfo_now.arena)/M;
#endif
}
