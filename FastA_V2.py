__author__ = 'Wombat'


class FastA:

    """
    V2.1  17 April 2020, Python 3.x conversion.
    This is a slightly smarter version of our previous FastA module. It is still not very careful about what
    properly constitutes a FastA sequence -- the big difference here is that we have made the class iterable by
    implementing the __iter__ and __next__ methods.
    In other words, I can now iterate over the fastA sequences present in some file, without having to first load
    each sequence into memory. The class takes an argument that is a filepath and each iteration will return the next
    (annotation, sequence) tuple in turn.  Probably this class could be improved by figuring out how to use the filehandle
    with a context manager -- right now the filehandle is used by multiple methods, so it's not clear to me how this
    might be accomplished.
    """

    translation_table = str.maketrans('', '', '0123456789*')

    def __init__(self, file_name):

        self.file_handle = open(file_name, 'r')
        self.annotation = ''

    def __iter__(self):

        """
        By defining an __iter__ method that just returns self, we are declaring this class to be an iterable
        Note that we can also do this using a generator / yield. But for now we will roll old-skool to see what goes on
        underneath the hood.
        """
        return self

    def __next__(self):

        """
        __next__ in Python 3.x is the key required method of an iterable class. It must simply return the
        next item in a sequence, or, if there are no more items, it should raise a StopIteration exception.
        Note that this implies that a normal loop in Python is actually ended, under-the-hood, by the error-handling
        mechanism! In the case of this next method, notice that the file handle is actually created in the __init__
        method and is held in an instance variable. That means that the current location of the read pointer in the
        file (indeed, the filehandle object itself) will remain persistent so long as the FastA object exists.
        """

        mode = 'scan'       # we always start off scanning for a new > character, indicating a new annotation
        sequence = ''       # We don't yet know what the sequence will be that we report back
        position = self.file_handle.tell()      # remember where we are in the file in case we need to rewind a line
        newline = self.file_handle.readline()

        while newline:
            newline = newline.strip()
            if mode == 'scan' and '>' in newline:   # We've discovered a new annotation. Grab it
                self.annotation = newline
                mode = 'collect'
                newline = self.file_handle.readline()
            elif mode == 'collect' and '>' in newline:  # We've finished collecting all the sequence
                self.file_handle.seek(position)         # pretend we didn't already read this line by rewinding
                return self.annotation[1:], sequence        # return the annotation and the sequence we discovered
            else:
                newline = newline.translate(self.translation_table)
                sequence += newline.upper()
                position = self.file_handle.tell()          # remember the current position in case we need to rewind
                newline = self.file_handle.readline()       # upon encountering the next annotation

        if self.annotation and sequence:                # mop up the last sequence
            return (self.annotation[1:], sequence)
        else:
            self.file_handle.seek(0)                # reset to the beginning of the file for another possible iteration
            raise StopIteration()

    def close(self):

        self.file_handle.close()


def main():

    """
    A quick demo of how a FastA object can be iterated across using a standard for loop.
    """

    sequence_file = FastA('6803Orfs.txt')

    for i, sequence in enumerate(sequence_file):

        print(i, sequence)

    sequence_file.close()


if __name__ == '__main__':

    main()
