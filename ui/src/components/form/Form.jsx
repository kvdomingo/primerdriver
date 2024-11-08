import {
  MDBContainer as Container,
  MDBIcon as Icon,
  MDBTypography as Typography,
} from "mdbreact";
import PropTypes from "prop-types";
import { Link } from "react-router-dom";

function Form(props) {
  return (
    <Container className="my-5">
      <Link className="btn btn-blue-grey btn-rounded mb-4" to="/" id="back">
        <Icon fas icon="arrow-left" className="mr-2" /> main menu
      </Link>
      <Typography
        tag="h2"
        variant="h2-responsive"
        className="d-md-inline ml-md-4 py-2 text-md-center"
      >
        {props.title}
      </Typography>
      <form
        id="form"
        className="form"
        onChange={props.handleValidate}
        onMouseUp={props.handleValidate}
        onKeyUp={props.handleValidate}
        onSubmit={props.handleSubmit}
      >
        {props.children}

        <div className="my-3 text-center">
          <input
            type="reset"
            id="reset"
            onClick={props.handleReset}
            className="btn btn-warning text-dark"
            value="Reset"
          />
          <input
            type="submit"
            id="submit"
            onClick={props.handleSubmit}
            className="btn btn-primary"
            value="Submit"
            disabled={!props.isValid}
          />
        </div>
      </form>
    </Container>
  );
}

Form.propTypes = {
  handleValidate: PropTypes.func,
  handleReset: PropTypes.func,
  handleSubmit: PropTypes.func,
  isValid: PropTypes.bool,
  title: PropTypes.string,
};

export { Form };
