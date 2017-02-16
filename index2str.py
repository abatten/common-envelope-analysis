def index2str(index_input, num_digits=4, append_character='0'):
        '''
        index2str converts an integer to a string to
         represent that integer.

        This is used so that the output files will be sorted in the
        correct order.

        i.e. 0000, 0001, 0002 etc rather than 1,10,11,12, etc

	index_input: Expects an integer.
	num_digits: Expects an integer.
	append_character: Expects a string.


        Examples:
        index2str(4) returns "0004"
        index2str(105) returns "0105"
	index2str(105, num_digits=6) returns "000105"
        index2str(40, append_character='t') returns "tt40"
	'''

        index_str = str(index_input)
        num_to_append = num_digits - len(index_str) # How many characters to add.
        new_str_index = [] # Stores the digits of the new index

        for i in range(num_to_append):
                new_str_index.append(append_character)

        new_str_index.append(index_str)
        index_str = ''.join(new_str_index)

        return(index_str)
