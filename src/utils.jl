using Compat.TypeUtils: typename

if :wrapper in fieldnames(TypeName)
    parameterless_type(T::Type) = typename(T).wrapper
else
    parameterless_type(T::Type) = typename(T).primary
end

parameterless_type(x) = parameterless_type(typeof(x))
