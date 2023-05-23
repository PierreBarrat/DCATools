"""
	const DEFAULT_AA_MAPPING

`"-ACDEFGHIKLMNPQRSTVWY"`.
"""
const DEFAULT_AA_MAPPING = "-ACDEFGHIKLMNPQRSTVWY"

"""
	DEFAULT_MAPPING(q::Int)

If `q==21`, return `DEFAULT_AA_MAPPING`. Otherwise, return some string.
"""
function DEFAULT_MAPPING(q::Int)
	return if q == 21
		DEFAULT_AA_MAPPING
	elseif q < 21
		DEFAULT_AA_MAPPING[1:q]
	else
		"A"^q
	end
end
