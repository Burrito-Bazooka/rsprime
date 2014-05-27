class Mapper():
    """
    alphabet is a string such as
    
    '0123456789abcdefghijkmnopqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ'

    equivs, if present is a list of tuples like this:

    [('O', '0'), ('I', '1'), ('l', '1')]
    """
    def __init__(self, alphabet, equivs=[]):
        self.table_c2n = {}
        self.table_n2c = {}
        for i in xrange(len(alphabet)):
            self.table_c2n[alphabet[i]] = i
            self.table_n2c[i] = alphabet[i]
        for tup in equivs:
            self.table_c2n[tup[0]] = self.table_c2n[tup[1]]

    def encode(self, data):
        return self._conv(self.table_n2c, data)

    def decode(self, data):
        return self._conv(self.table_c2n, data)

    def pad(self, s, w):
        return s.rjust(w, self.table_n2c[0])

    def strip(self, s):
        return s.lstrip(self.table_n2c[0])

    def _conv(self, table, data):
        if isinstance(data, int):
            if data in table:
                return table[data]
            else:
                raise ValueError("No entry found in conversion table for '%r'!" % data)
        if isinstance(data, str):
            if len(data) == 1:
                if data in table:
                    return table[data]
                else:
                    return -1
            else:
                return self._conv(table, list(data))
        if isinstance(data, (list, tuple)):
            r = map(lambda x: self._conv(table, x), data)
            if isinstance(r[0], str):
                return ''.join(r)
            else:
                return r

        raise ValueError("Unsupported type '%s'!" % type(data)) 

# vim: sw=4 ts=4 et ai si bg=dark
