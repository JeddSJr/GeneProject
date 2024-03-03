import { Component } from '@angular/core';
import {FormBuilder, Validators} from '@angular/forms';
import {DropdownModule} from 'primeng/dropdown'
@Component({
  selector: 'app-user-form',
  templateUrl: './user-form.component.html',
  styleUrl: './user-form.component.css'
})

export class UserFormComponent {
  ethnicityOptions = ['Asian', 'African','Caucasian']
  userForm = this.fb.group({ 
    firstName: ['', Validators.required],
    lastName: ['', Validators.required],
    ethnicity: ['', Validators.required],
    height : ['', Validators.required],
    weight : ['', Validators.required]
  });
  constructor(public fb: FormBuilder){
  
  }
  submit(){

  }
}
