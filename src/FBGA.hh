/*
(***********************************************************************************)
(*                                                                                 *)
(* The FBGA project                                                               *)
(*                                                                                 *)
(* Copyright (c) 2025, Mattia Piazza                                               *)
(*                                                                                 *)
(*    Mattia Piazza                                                                *)
(*    Department of Industrial Engineering                                         *)
(*    University of Trento                                                         *)
(*    e-mail: mattia.piazza@unitn.it                                               *)
(*                                                                                 *)
(***********************************************************************************)
*/

///
/// file: FBGA.hh
///

#pragma once

#ifndef INCLUDE_FBGA
#define INCLUDE_FBGA


// Print FBGA errors
#ifndef FBGA_ERROR
#define FBGA_ERROR(MSG)                                                                           \
  {                                                                                                \
    throw std::runtime_error(std::to_string(MSG));                                                 \
  }
#endif

// Check for FBGA errors
#ifndef FBGA_ASSERT
#define FBGA_ASSERT(COND, MSG)                                                                    \
  if (!(COND))                                                                                     \
  FBGA_ERROR(MSG)
#endif


#endif

///
/// eof: FBGA.hh
///
