import { Link } from "@tanstack/react-router";
import { MDBContainer as Container, MDBIcon as Icon } from "mdbreact";
import type { ReactNode } from "react";
import { usePrimerDriverStore } from "../../store";

interface FormProps {
  onSubmit: () => void;
  isValid: boolean;
  title: string;
  children: ReactNode;
}

function Form(props: FormProps) {
  const resetForm = usePrimerDriverStore(state => state.resetForm);

  return (
    <Container className="my-5">
      <Link className="btn btn-blue-grey btn-rounded mb-4" to="/" id="back">
        <Icon fas icon="arrow-left" className="mr-2" /> main menu
      </Link>
      <h2 className="h2-responsive d-md-inline ml-md-4 py-2 text-md-center">
        {props.title}
      </h2>
      <form
        id="form"
        className="form"
        onSubmit={e => {
          e.preventDefault();
          props.onSubmit();
        }}
      >
        {props.children}

        <div className="my-3 text-center">
          <input
            type="reset"
            id="reset"
            onClick={resetForm}
            className="btn btn-warning text-dark"
            value="Reset"
          />
          <input
            type="submit"
            id="submit"
            className="btn btn-primary"
            value="Submit"
            disabled={!props.isValid}
          />
        </div>
      </form>
    </Container>
  );
}

export { Form };
