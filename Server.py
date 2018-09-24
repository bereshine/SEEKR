#  all the following function will be called by parse_mail_from_cmd
#  I tiried to report the MAIL error as deep as possible.
#  all the following function besides parse_mail_from_cmd will either return True,current(position of string)
#  or return False,0.5 (0.5 is just a flag to make the mistake being able to be printed just once by the
#  deepest function that found the mistake)


#HW2 : the state machine automata is as following:
# n--(conneciton)--o--(helo)---p--(mail from)--q--(rcpt)---r--(data/rcpt)--s(content)--(s)

import socket
def parse_quit(string, current):

    if len(string)>=4 and string[current:current+4]== 'QUIT':
        current=current+4
        bolean,current = parse_null_space(string,current)
        if current == len(string):
            return True
        else:
            return False

    else:
        return False

def parse_helo(string, current):

    if len(string)>=4 and string[current:current+4] == "HELO":
        current=current+4
        bolean, current = parse_domain(string,current)
        if bolean == True:
            return 'p'
        else:
            return 'n'
    else:
        return 'n'


def parse_content(string, current):

        if string == '.':
        # print "250 OK"
            return 'p'
        else:

            return 's'


def parse_data(string, current):

    if len(string)>=4 and string[current:current+4] == "DATA":
        current=current+4
        bolean, current = parse_null_space(string, current)
        if current<len(string):
            if current==4:

                return "p500"
            else:

                return "p501"
        else:
            return "s"
    else:

        return "p500"


def parse_rcpt_to_cmd(string, current):
    if string[0:4] == "RCPT":
        bolean, current = parse_white_space(string, 4)
        if bolean == True and string[current:current + 3] == "TO:":
            current = current + 3;
            bolean2, current = parse_null_space(string, current)
            if current < len(string) and bolean2 == True:
                bolean3, current = parse_reverse_path(string, current)
                if bolean3 == True:
                    bolean4, current = parse_null_space(string, current)
                    if string[current:len(string)] == "":
                        # print "250 OK"
                        return "r"
                    else:
                        # print "501 Syntax error in parameters or arguments"
                        return "p501"
                elif current == 0.5:
                    return "p501"
                else:
                    # print "501 Syntax error in parameters or arguments"
                    return "p501"
            else:
                # print "501 Syntax error in parameters or arguments"
                return "p501"
        elif current == 0.5:
            return "p500"
        else:
            # print "500 Syntax error: command unrecognized"
            return "p500"
    else:
        # print "500 Syntax error: command unrecognized"
        return "p500"


def parse_mail_from_cmd(string, current):
    if string[0:4] == "MAIL":
        bolean, current = parse_white_space(string, 4)
        if bolean == True and string[current:current + 5] == "FROM:":
            current = current + 5;
            bolean2, current = parse_null_space(string, current)
            if current < len(string) and bolean2 == True:
                bolean3, current = parse_reverse_path(string, current)
                if bolean3 == True:
                    bolean4, current = parse_null_space(string, current)
                    if string[current:len(string)] == "":
                        # print "250 OK"

                        return "q"
                    else:
                        # print "501 Syntax error in parameters or arguments"
                        return "p501"
                elif current == 0.5:
                    return "p501"
                else:
                    # print "501 Syntax error in parameters or arguments"
                    return "p501"
            else:
                # print "501 Syntax error in parameters or arguments"
                return "p501"
        elif current == 0.5:
            return "p500"
        else:
            # print "500 Syntax error: command unrecognized"
            return "p500"
    else:
        # print "500 Syntax error: command unrecognized"
        return "p500"


def parse_white_space(string, current):
    if string[current] != " " and string[current] != "\t":
        # print "500 Syntax error: command unrecognized"
        return False, 0.5
    else:
        while string[current] == " " or string[current] == "\t":
            current = current + 1
        return True, current


def parse_null_space(string, current):
    while current < len(string) and (string[current] == " " or string[current] == "\t"):
        current = current + 1
    return True, current


def parse_reverse_path(string, current):
    bolean, current = parse_path(string, current)
    return bolean, current


def parse_path(string, current):
    if current < len(string) and string[current] == "<":
        current = current + 1
        bolean, current = parse_mailbox(string, current)
        if current < len(string) and bolean == True and string[current] == ">":
            current = current + 1
            return True, current
        elif current == 0.5:
            return False, 0.5
        else:
            # print "501 Syntax error in parameters or arguments"
            return False, 0.5
    else:
        # print "501 Syntax error in parameters or arguments"
        return False, 0.5


def parse_mailbox(string, current):
    if current < len(string):
        bolean, current = parse_local_part(string, current)
        if bolean == True and string[current] == "@":
            bolean2, current = parse_domain(string, current + 1)
            if bolean2 == True:
                return True, current
            elif current == 0.5:
                return False, 0.5
            else:
                # print "501 Syntax error in parameters or arguments"
                return 0
        elif current == 0.5:
            return False, 0.5
        else:
            # print "501 Syntax error in parameters or arguments"
            return False, 0.5
    else:
        # print "501 Syntax error in parameters or arguments"
        return False, 0.5


def parse_local_part(string, current):
    if current < len(string):
        bolean, current = parse_string(string, current)
    return bolean, current


def parse_string(string, current):
    special = '<>()[]\.,:;@" \t'
    if current + 1 < len(string):
        bolean = parse_char(string[current])
        if bolean == 0.5:
            return False, 0.5
        if bolean == True and (0 <= ord(string[current + 1]) <= 127 and string[current + 1] not in special) != True:
            current = current + 1
            return True, current
        elif bolean == True and (0 <= ord(string[current + 1]) <= 127 and string[current + 1] not in special) == True:
            bolean2, current = parse_string(string, current + 1)
            return bolean2, current
        else:
            # print "501 Syntax error in parameters or arguments"
            return False, 0.5
    else:
        # print "501 Syntax error in parameters or arguments"
        return False, 0.5


def parse_char(string):
    special = '<>()[]\.,:;@" \t'
    if 0 <= ord(string) <= 127 and string not in special:
        #             ord(string) != 60 and  ord(string) != 62 and  ord(string) !=40 \
        # and ord(string)!=41 and ord(string)!=91 and ord(string)!=93 and ord(string)!=92 and ord(string) !=46 \
        # and  ord(string)!= 44 and ord(string)!=58 and ord(string)!=59 and  ord(string)!=64 and ord(string)!=34:
        # #current=current+1
        return True
    else:
        # print "501 Syntax error in parameters or arguments"
        return 0.5


def parse_domain(string, current):
    bolean, current = parse_element(string, current)
    if current == 0.5:
        return False, 0.5
    if bolean == True and string[current] == ".":
        bolean2, current = parse_domain(string, current + 1)
        return bolean2, current
    elif bolean == True:
        return bolean, current
    else:
        # print"501 Syntax error in parameters or arguments"
        return False, 0.5


def parse_element(string, current):
    bolean, current = parse_name(string, current)
    return bolean, current


def parse_name(string, current):
    bolean = parse_alpha(string[current])
    if bolean == True and (65 <= ord(string[current + 1]) <= 90 or 97 <= ord(string[current + 1]) <= 122 or 48 <= ord(
            string[current + 1]) <= 57):
        bolean2, current = parse_let_dig_str(string, current + 1)
        return True, current
    elif bolean == 0.5:
        return False, 0.5
    else:
        # print "501 Syntax error in parameters or arguments"
        return False, 0.5


def parse_let_dig_str(string, current):
    if (65 <= ord(string[current]) <= 90 or 97 <= ord(string[current]) <= 122 or 48 <= ord(string[current]) <= 57):
        if current < len(string) - 1:
            current = current + 1

            bolean, current = parse_let_dig_str(string, current)
            return True, current
        else:
            return True, current
    elif (65 <= ord(string[current]) <= 90 or 97 <= ord(string[current]) <= 122 or 48 <= ord(
            string[current]) <= 57) != True:

        return True, current
    else:
        current = current + 1
        return True, current


def parse_let_dig(string):
    if parse_alpha(string) or parse_digit(string):
        return True
    else:
        # print "501 Syntax error in parameters or arguments"
        return 0


def parse_alpha(string):
    if 65 <= ord(string) <= 90 or 97 <= ord(string) <= 122:
        return True
    else:

        # print "501 Syntax error in parameters or arguments"
        return 0.5


def parse_digit(string):
    if 48 <= ord(string) <= 57:
        return True
    else:
        # print "501 Syntax error in parameters or arguments"
        return 0

# initialization
import sys
serverPort=int(sys.argv[1])


# serverPort=8000+4563

serverSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
serverSocket.bind(('', serverPort))
serverSocket.listen(1)

state = 'n'
rcpt_list = []
text = []
while True:
    try:
        # save keyboard input

        # input = raw_input()

        if state == 'n':
            connectionSocket, addr=serverSocket.accept()
            # input = connectionSocket.recv(1024).decode()
            # conn_res='220'+socket.gethostname()
            conn_res='220'+socket.gethostname()

            connectionSocket.send(conn_res.encode())
            state = 'o'
        elif state == 'o':
            input = connectionSocket.recv(1024).decode()

            helo_state = parse_helo(input,0)
            if helo_state == 'p':
                state = 'p'
                helo_res='250'+input+'pleased to meet you.'
                connectionSocket.send(helo_res.encode())
            else:
                str = 'HELO error'
                print str
                connectionSocket.send(str.encode())
                state = 'n'

        elif state == 'p':
            input = connectionSocket.recv(1024).decode()

            mail_state = parse_mail_from_cmd(input, 0)
            rcpt_state = parse_rcpt_to_cmd(input, 0)
            data_state = parse_data(input, 0)
            if parse_quit(string,current)==False:
                state = 'n'
            if mail_state == "q":
                str = '250 OK'

                connectionSocket.send(str.encode())
                state = 'q'

                start_mail = input.find('<')
                end_mail = input.find('>')+1
                mail = input[start_mail:end_mail]
                # print mail
            elif mail_state == 'p500':
                if rcpt_state == 'p501'or rcpt_state=='r' or data_state == 'p501' or data_state == 's':


                    str = '503 Bad sequence of commands'
                    print str
                    connectionSocket.send(str.encode())
                    state = 'n'
                elif rcpt_state == 'p500' or data_state == 'p500':



                    str ='500 Syntax error: command unrecognized'
                    print str
                    connectionSocket.send(str.encode())
                    state = 'n'

            else:
                str = '501 Syntax error in parameters or arguments'

                print str
                connectionSocket.send(str.encode())
                state = 'n'

        elif state == 'q':
            input = connectionSocket.recv(1024).decode()

            mail_state = parse_mail_from_cmd(input, 0)
            rcpt_state = parse_rcpt_to_cmd(input, 0)
            data_state = parse_data(input, 0)

            if rcpt_state == "r":

                str = '250 OK'
                state = 'r'

                start_rcpt = input.find('<')
                end_rcpt = input.find('>')+1
                rcpt = input[start_rcpt:end_rcpt]
                rcpt_list.append(rcpt)

            elif rcpt_state == 'p500':
                if mail_state == 'p501' or mail_state=='q' or data_state == 'p501' or data_state == 's':

                    str= '503 Bad sequence of commands'
                    print str
                    connectionSocket.send(str.encode())
                    state = 'n'
                    rcpt_list = []
                elif mail_state == 'p500' or data_state == 'p500':
                    str = '500 Syntax error: command unrecognized'
                    print str
                    connectionSocket.send(str.encode())
                    state = 'n'
                    rcpt_list = []

            else:
                str ='501 Syntax error in parameters or arguments'
                print str
                connectionSocket.send(str.encode())
                state = 'n'
                rcpt_list = []
        elif state == 'r':
            input = connectionSocket.recv(1024).decode()

            mail_state = parse_mail_from_cmd(input, 0)
            rcpt_state = parse_rcpt_to_cmd(input, 0)
            data_state = parse_data(input, 0)

            if rcpt_state == 'r':
                str = '250 OK'
                connectionSocket.send(str.encode())
                state = 'r'

                start_rcpt = input.find('<')
                end_rcpt = input.find('>')+1
                rcpt = input[start_rcpt:end_rcpt]

                rcpt_list.append(rcpt)
                # print rcpt_list
            elif data_state == 's':

                str = '354 Start mail input; end with <CRLF>.<CRLF>'
                connectionSocket.send(str.encode())
                state = 's'
            elif rcpt_state == 'p500' and data_state == 'p500':
                if mail_state == 'p501' or mail_state == 'q':



                    str = '503 Bad sequence of commands'
                    print str
                    connectionSocket.send(str.encode())

                    state = 'n'
                    rcpt_list = []
                elif mail_state == 'p500':

                    str = '500 Syntax error: command unrecognized'
                    print str
                    connectionSocket.send(str.encode())
                    state = 'n'
                    rcpt_list = []

            else:

                str = '501 Syntax error in parameters or arguments'
                print str
                connectionSocket.send(str.encode())
                state = 'n'
                rcpt_list = []

        elif state == 's':
            input = connectionSocket.recv(1024).decode()


            state = parse_content(input, 0)
            if state == 'p':
                str = '250 OK'
                connectionSocket.send(str.encode())
                for forward_path in rcpt_list:
                    left_domain = forward_path.find('@')+1

                    file = open("forward/" + forward_path[left_domain:len(forward_path)], "a+",0)
                    file.write("From: " + mail + "\n")
                    for to in rcpt_list:
                        file.write("To: " + to + "\n")
                    for line in text:
                        file.write(line + "\n")
                    file.close

                rcpt_list = []
                text = []

            else:

                text.append(input)
    except (EOFError):
        break  # end of file reached




