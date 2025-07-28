
def setup_env():





    return envs

def setup():
    whole_dic={}
    dic=setup_env()
    whole_dic.update({'global':dic})
    
    env_file='environment_for_server_config.yaml'
    with open(env_file, 'w') as wr:
        yaml.safe_dump(whole_dic, wr)
        string=yaml.safe_dump(whole_dic)



if __name__=='__main__':
    setup()
