I provide the base for preconditioners. All preconditioners should implement `applyDirect:` and `applyInverse:` methods of multiplying
the preconditioner by a vector.