from chrmapper import ChrMapper

if __name__ == "__main__":
    alphabet = '0123456789abcdefghijkmnopqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ'
    equivs = [('O', '0'), ('I', '1'), ('l', '1')]

    mapper = ChrMapper(alphabet, equivs)
    dec = mapper.decode('AllYourBaseAreBelongToUs')
    enc = mapper.encode(dec)
    print dec
    print enc
    dec = mapper.decode('0123456789abcdeghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    enc = mapper.encode(dec)
    print dec
    print enc
    print mapper.encode(range(len(alphabet)))
