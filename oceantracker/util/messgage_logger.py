from os import path, remove
import traceback

class GracefulError(Exception):
    def __init__(self, message='-no error message given',hint=None):
        # Call the base class constructor with the parameters it needs
        msg= 'Error >> ' + message + '\n hint= ' + hint if hint is not None else ' Look at messages above or in .err file'
        super(GracefulError, self).__init__(msg)

def msg_str(msg,tabs=0):
    tab = '  '
    m = ''
    for n in range(tabs): m += tab
    m +=  msg
    return m

class MessageLogger(object):
    def __init__(self,screen_tag,max_warnings = 50):
        self.screen_tag = screen_tag
        self.max_warnings = max_warnings
        self.fatal_error_count = 0
        self.warnings_and_errors=[]
        self.log_file = None
        self.error_warning_count = 0

    def set_up_files(self,run_output_dir,output_file_base):

        # log file
        log_file_name = output_file_base + '_log.txt'
        self.log_file_name= path.join(run_output_dir, log_file_name)
        self.log_file= open(self.log_file_name, 'w')

        # kill any old error file
        error_file_name = output_file_base+ '.err'
        self.error_file_name = path.join(run_output_dir, error_file_name)
        if path.isfile(self.error_file_name ):
            remove(self.error_file_name)

        return  log_file_name, error_file_name
    #todo add abilty to return excecption/traceback?
    def msg(self, msg_text, warning=False, note=False,
            hint=None, tag=None, tabs=0, crumbs=None,
            fatal_error=False,exit_now=False, exception = None, traceback_str=None):

        if exception is not None:
            fatal_error = True
            exit_now = True

        if fatal_error: self.fatal_error_count +=1

        m = ['']

        # first line of message
        if fatal_error:
            m[0] += msg_str( '>>> Error: ', tabs)
            self.error_warning_count += 1

        elif warning:
            m[0] += msg_str('>>> Warning: ' , tabs)
            self.error_warning_count += 1
        elif note:
            m[0] += msg_str('>>> Note: ', tabs)
        else:
            m[0] += msg_str('', tabs)

        if exception is not None:
            m[0] += msg_str('exception >>: ' + str(exception), tabs+2)

        if traceback_str is not None:
            m[0] += msg_str('traceback >>: ' + str(traceback_str), tabs + 2)

        m[0] +=  msg_text

        # first line complete

        if crumbs is not None:
            m.append(msg_str('in: ' + crumbs, tabs + 3))
        if hint is not None:
            m.append(msg_str('Hint: ' + hint, tabs + 3))

        # write message lines
        for l in m:
            print(self.screen_tag + ' ' + l)
            if self.log_file is not None:
                self.log_file.write(l + '\n')

            # keeplist ond warnings errors etc to print at end
            if fatal_error or warning or note:
                if self.error_warning_count <= self.max_warnings:
                    self.warnings_and_errors.append(l)

        # mange whether to exit now or not:
        if fatal_error: self.fatal_error_count += 1

        # todo add traceback to message?
        if exit_now:
            raise GracefulError('Fatal error cannot continue')

    def has_fatal_errors(self): return  self.fatal_error_count > 0

    def exit_if_prior_errors(self,msg=None):
        if self.has_fatal_errors():
            raise GracefulError('Fatal error cannot continue >>> ' +msg if msg is not None else '', hint='Check above or run.err file for errors')

    def insert_screen_line(self):
        self.msg('--------------------------------------------------------------------------')

    def progress_marker(self, msg, tabs=0):
        self.msg('- ' + msg, tabs=tabs + 1)

    def show_all_warnings_and_errors(self,error=None):
        for l in self.warnings_and_errors:
            print(self.screen_tag + ' ' + l)
            if self.log_file is not None:
                self.log_file.write(l + '\n')
        if error is not None:
            self.msg(str(error))
            self.msg(traceback.format_exc())



    def write_error_log_file(self, e=None):

        if self.log_file is None : return

        with open(path.normpath(self.error_file_name),'w') as f:
            f.write('_____ Known warnings and Errors ________________________________\n')
            for l in self.warnings_and_errors:
                f.write(l + '\n')


            f.write('________Trace back_____________________________\n')

            if e is not None:
                f.write(str(e))
                self.msg(str(e))
                s = traceback.format_exc()
                f.write(s)
                self.msg(s)

    def close(self):
        if self.log_file is not None: self.log_file.close()