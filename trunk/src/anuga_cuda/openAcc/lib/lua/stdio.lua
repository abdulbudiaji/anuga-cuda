
-- A simple implementation of the most common functions from <stdio.h>.
-- using the standard io package from lua.
--
-- The following functions and variables are provided 
--
-- 
--  function fprintf(file, fmt, ...)         Format to a 'file'  
--  function printf(fmt,...)                 alias for fprintf(stdout, fmt, ...)
--  function sprintf(fmt,...)                Format to a string  
--  function file = fopen(filename, mode)    Open a 'file'. 
--                                           'mode' can be one of "r", "w", "a", "r+", "w+" or "a+"
--                                           with optional "b" at the end for binary files (on 
--                                           some systems)
--  function fclose(file)                    Close a file 
--  function fflush(file)                    Flush a file to disk
--  function fgetc(file)                     Read a single character from 'file' 
--  function fgets(file)                     Read a single line from 'file' (api is different from <stdio.h>)
--  
--  function fseek(file,whence,offset)       Set the position in a file.
--                                           Argument whence can be one of 
--                                            SEEK_SET or "set"
--                                            SEEK_CUR or "cur"
--                                            SEEK_END or "end"
--              
--  function rewind(file)                    Equivalent to fseek("set",0) 
--  
--  stdin                           The standard input file 
--  stderr                          The error output file 
--  stdout                          The error output file 
--
-- The files used by those function are compatible with those 
-- provided by the default lua packages 'io' and 'file' 
--
-- The following functions are missing 
--   scanf() ...   :  use the standard lua functions such as str:match() instead
--   feof()        :  user should check that the result of fgetc() or fgets() is not nul 
--
-- 

function printf(fmt,...)
   return io.write(fmt:format(...))
end 

function sprintf(fmt,...)
   return fmt:format(...) 
end 

function fopen(filename, mode)
   return io.open(filename,mode)
end 

function fclose(file)
   return file:close() 
end 

function fflush(file)
   return file:flush() 
end 

function fgetc(file)
   return file:read(1) 
end 

--
-- Read a single line 
--
-- Warning: the api is not similar to the C version
-- 
function fgets(file)
   return file:read() 
end 

function fprintf(file , fmt, ...)
   return file:write( fmt:format(...) )
end 


stdin  = io.stdin 
stderr = io.stderr 
stdout = io.stdout
 
SEEK_SET="set"
SEEK_CUR="cur"
SEEK_END="end"
