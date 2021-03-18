def standard_encode(i,N):
    format_style = '0' + str(N) + 'b'
    return format(i, format_style)


def gray_code(i,N):
    gray=i^(i>>1)
    return "{0:0{1}b}".format(gray,N)