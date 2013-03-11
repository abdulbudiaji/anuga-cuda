HEADER = (
#include \<math.h>
#include \<stdio.h>
#include \<stdlib.h>

#include \<hmpprt/HostTypes.h>
#include \<hmpprt/Context.h>
#include \<hmpprt/Grouplet.h>
#include \<hmpprt/DeviceManager.h>
#include \<hmpperr/hmpperr.h>

#ifdef _WIN32
#  define CDLT_API __declspec(dllexport)
#else /* ! _WIN32 */
#  define CDLT_API
#endif /* _WIN32 */

<forward_declarations>

<dwarf_xml>

<host_type_definitions>
)

SIMD_LENGTH = #define HMPPCG_SIMD_LENGTH <value>\n\n

FUNCTION = <linkage end=" "><type> <name>(<params separator=", ">)\n

AFTER_FUNCTION = \n\n\n

BEGIN_BLOCK = <before>{\n <prologue end="\n">

END_BLOCK = <epilogue end="\n">}\n<after>

STRUCT_FORWARD = struct <name>;\n

STRUCT = (struct <name>
{
  <fields separator=";\n" end=";">
};\n)

UNION_FORWARD = union <name>;\n

UNION = (union <name>
{
  <fields separator=";\n" end=";">
};\n)

IF = if (<cond>)\n

ELSE = else\n

FOR = for (<init separator=", "> ; <cond> ; <cycle separator=", ">)\n

WHILE = while (<cond>)\n

POSTWHILE = while (<cond>);\n

DO = do\n

BREAK = break;\n

CONTINUE = continue;\n

RETURN = return <expr>;\n

GOTO = goto <name>;\n

LABEL = <name>:;\n

DECLARE = <declare>;\n

INIT_SIGNATURE = hmpprt::Context::getInstance()->getGrouplet()->addSignature("<name>", "<value>");\n

SET_TARGET = hmpprt::Context::getInstance()->getGrouplet()->setTarget(hmpprt::<value>);\n

HMPPRT_INIT = (
extern "C" CDLT_API void * hmpprt_init()
{
  try
  {
    <first_time_content>
    <every_time_content>
  }
  catch (hmpperr::Error & e)
  {
    return e.clone();
  }
  catch (...)
  {
    fprintf(stderr, "Unexpected error in hmpprt_init()\\n");
    abort();
  }
  return 0;
}
)

HMPPRT_FINI = (
extern "C" CDLT_API void * hmpprt_fini()
{
  try
  {
    <content>
  }
  catch (hmpperr::Error & e)
  {
    return e.clone();
  }
  catch (...)
  {
    fprintf(stderr, "Unexpected error in hmpprt_fini()\\n");
    abort();
  }
  return 0;
}
)

DWARFXML = (
static const char __hmpp_dwarf_xml[] = 
<cstring>;
)

FOOTER = (
// footer
)
