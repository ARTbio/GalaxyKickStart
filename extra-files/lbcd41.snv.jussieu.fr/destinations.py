from galaxy.jobs.mapper import JobMappingException

def cluster_heavy(user, tool):
    # user is a galaxy.model.User object or None
    # tool is a galaxy.tools.Tool object
    if user is None:
        raise JobMappingException('You must login to use this tool!')
    else:
        required_role = "IBPS"
        final_destination = "cluster_heavy_ibps"
    # Check that the required_role is in the set of role names associated with the user
    user_roles = user.all_roles() # a list of galaxy.model.Role objects
    user_in_role = required_role in [role.name for role in user_roles]
    if not user_in_role:
        return "cluster_heavy"
    else:
        return final_destination

def cluster(user, tool):
    # user is a galaxy.model.User object or None
    # tool is a galaxy.tools.Tool object
    if user is None:
        return "cluster"
    else:
        required_role = "IBPS"
        final_destination = "cluster_ibps"
    # Check that the required_role is in the set of role names associated with the user
    user_roles = user.all_roles() # a list of galaxy.model.Role objects
    user_in_role = required_role in [role.name for role in user_roles]
    if not user_in_role:
        return "cluster"
    else:
        return final_destination


