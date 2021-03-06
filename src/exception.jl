#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

import Base: showerror

export
    DescriptorError,
    OperationError


struct DescriptorError <: Exception
    msg::String
end

showerror(io::IO, e::DescriptorError) = print(io, e.msg)


struct OperationError <: Exception
    msg::String
end

showerror(io::IO, e::OperationError) = print(io, e.msg)
